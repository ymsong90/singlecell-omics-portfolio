################################################################################
# Serial Slide Image Registration - Napari Interactive Method
# reg_serial.py
#
# Description:
#   Interactive image registration between IMC (ROI) and Fontana-Masson
#   stained tissue slides (Slide) using napari GUI. Supports manual transform
#   with real-time visualization and pixel-wise remapping for high precision.
#
# Features:
#   - Manual transform control (translate, rotate, scale)
#   - Real-time ROI boundary visualization
#   - Numeric rotation input for precise angles
#   - Pixel-wise remapping using cv2.remap
#   - DPI-preserved TIFF output
#   - Batch processing support
#
# Input:
#   - ROI files: *_XXX_YYY.tif (e.g., NHS_001_001.tif)
#   - Slide files: *.tif (e.g., NHS_001.tif)
#
# Output:
#   - roi.tif: Original ROI image
#   - registered.tif: Pixel-remapped aligned slide
#   - overlay.tif: RGB overlay (R=ROI, G=registered)
#
# Usage:
#   Single pair:
#     python reg_serial.py single --roi roi.tif --slide slide.tif --out results
#   
#   Batch processing:
#     python reg_serial.py batch --dir input_folder --out results
#
# Author: YMS
# Organization: Seoul National University Hospital, Mokam Research Institute
# Date: 2025-05-01
################################################################################

import sys
import argparse
import re
from pathlib import Path

import numpy as np
import cv2
from skimage import io
from PIL import Image

################################################################################
# 1. File Discovery and Pairing Functions
################################################################################

def load_images(folder):
   
    from pathlib import Path
    import re

    folder = Path(folder)
    roi_files = []
    slide_files = []
    
    for f in sorted(folder.glob('*.tif')):
        # ROI files must have _XXX_YYY.tif pattern
        # Example: NHS_001_001.tif, sample_003_002.tif
        if re.search(r'_[0-9]{3}_[0-9]{3}\.tif$', f.name):
            roi_files.append(f)
        else:
            slide_files.append(f)
    
    return roi_files, slide_files


def strip_suffix(p: Path):

    return re.sub(r'_[0-9]{3}$', '', p.stem)


def match_pairs(folder: Path):

    roi_files, slide_files = load_images(folder)
    pairs = {}
    
    for roi in roi_files:
        base = strip_suffix(roi)
        
        # Try exact match
        exact = folder / f"{base}.tif"
        if exact in slide_files:
            pairs[roi] = exact
            continue
        
        # Try numeric prefix match
        m = re.match(r'^(\d+)', base)
        if not m:
            continue
        
        sid = m.group(1)
        candidate = folder / f"{sid}.tif"
        
        if candidate in slide_files:
            pairs[roi] = candidate
            continue
        
        # Try closest match with same numeric prefix
        cands = [s for s in slide_files if s.stem.startswith(f"{sid}_")]
        if cands:
            pairs[roi] = min(cands, key=lambda s: abs(len(s.stem) - len(sid)))
    
    # Print matching results
    print(f"▶ 매핑 {len(pairs)}/{len(roi_files)}")
    for r, s in pairs.items():
        print(f"  {r.name} -> {s.name}")
    
    # Report unmatched ROI files
    missing = set(roi_files) - set(pairs.keys())
    if missing:
        print("⚠️ 매칭 실패한 ROI:")
        for r in missing:
            print("   ", r.name)
    
    return pairs


def load_rgb(path):

    # Try OpenCV first (faster for most formats)
    img = cv2.imread(str(path))
    if img is not None:
        return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    
    # Fallback to skimage for other formats
    img = io.imread(path)
    
    if img.ndim == 2:
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB)
    elif img.shape[2] == 4:
        img = img[..., :3]
    
    return img.astype(np.uint8)

################################################################################
# 2. Napari Interactive Registration
################################################################################

def align_and_crop(roi_img, slide_img, slide_path, output_dir):

    try:
        import napari
        from qtpy.QtWidgets import (
            QWidget, QVBoxLayout, QPushButton,
            QLabel, QSlider, QDoubleSpinBox
        )
        from qtpy.QtCore import Qt
    except ImportError:
        print("❌ Missing dependencies. Install with:")
        print("   pip install napari[all] qtpy")
        return

    h, w = roi_img.shape[:2]
    result = {'saved': False, 'image': None}

    # Create napari viewer
    viewer = napari.Viewer(title="ROI 정합 & 숫자 회전 크롭")

    # ========================================
    # Layer 1: Fixed ROI with boundary outline
    # ========================================
    viewer.add_image(
        roi_img, 
        name="ROI (fixed)", 
        opacity=1.0, 
        colormap='gray'
    )
    
    # Add ROI boundary as red polygon
    # NOTE: Coordinates are (row, col) format
    outline = [(0, 0), (h-1, 0), (h-1, w-1), (0, w-1)]
    viewer.add_shapes(
        [outline],
        shape_type="polygon",
        edge_color="red",
        face_color="transparent",
        edge_width=2,
        name="ROI_outline"
    )

    # ========================================
    # Layer 2: Transformable slide
    # ========================================
    slide_layer = viewer.add_image(
        slide_img, 
        name="Slide (transform)",
        opacity=0.7, 
        blending="additive"
    )
    
    # IMPORTANT: Set to transform mode for manual alignment
    viewer.layers.selection.active = slide_layer
    viewer.layers.selection.active.mode = 'transform'

    # ========================================
    # Control Panel: Opacity and Rotation
    # ========================================
    ctrl = QWidget()
    layout = QVBoxLayout(ctrl)

    # Opacity slider
    layout.addWidget(QLabel("Slide 투명도"))
    opa_slider = QSlider(Qt.Horizontal)
    opa_slider.setRange(0, 100)
    opa_slider.setValue(int(slide_layer.opacity * 100))
    opa_slider.valueChanged.connect(
        lambda v: setattr(slide_layer, 'opacity', v / 100)
    )
    layout.addWidget(opa_slider)

    # Numeric rotation input
    layout.addWidget(QLabel("회전 각도 (°)"))
    angle_spin = QDoubleSpinBox()
    angle_spin.setRange(-180.0, 180.0)
    angle_spin.setSingleStep(0.1)
    angle_spin.setValue(0.0)
    layout.addWidget(angle_spin)

    def apply_numeric_rotation(angle_deg):
        """
        Apply rotation to slide layer using numeric input.
        
        Rotation is applied around the image center (w/2, h/2) in world
        coordinates. The rotation matrix is combined with existing transform.
        
        Args:
            angle_deg: Rotation angle in degrees (counterclockwise)
        """
        tf = slide_layer.affine
        M = getattr(tf, 'affine_matrix', np.array(tf))
        
        # Rotation center (world coordinates)
        cx, cy = w / 2, h / 2
        
        # Create rotation matrix
        R = cv2.getRotationMatrix2D((cx, cy), angle_deg, 1.0)
        M2 = np.eye(3)
        M2[:2, :] = R
        
        # Combine with existing transform
        tf.affine_matrix = M2 @ M

    angle_spin.valueChanged.connect(apply_numeric_rotation)

    # ========================================
    # Save Button and Keyboard Shortcut
    # ========================================
    def on_save():
        """
        Save registered image using pixel-wise remapping.
        
        Algorithm:
            1. Generate pixel grid for ROI dimensions (h x w)
            2. Map each ROI pixel to slide coordinates using napari transform
            3. Apply cv2.remap for high-precision resampling
            4. Fill empty regions with edge pixels (BORDER_REPLICATE)
        
        NOTE: This approach preserves sub-pixel accuracy from napari transform
        """
        # Generate ROI pixel grid
        yy, xx = np.mgrid[0:h, 0:w]
        
        # Map world coordinates to data coordinates
        # IMPORTANT: Use slide_layer.world_to_data for accurate mapping
        mapped = np.array([
            slide_layer.world_to_data((float(r), float(c)))
            for r, c in zip(yy.ravel(), xx.ravel())
        ], dtype=float)
        
        map_y = mapped[:, 0].reshape(h, w).astype(np.float32)
        map_x = mapped[:, 1].reshape(h, w).astype(np.float32)

        # Apply remapping with linear interpolation
        cropped = cv2.remap(
            slide_img, 
            map_x, 
            map_y,
            interpolation=cv2.INTER_LINEAR,
            borderMode=cv2.BORDER_REPLICATE  # Fill gaps with edge pixels
        )

        result['saved'] = True
        result['image'] = cropped
        viewer.close()

    btn = QPushButton("저장 & 종료 (S)")
    btn.clicked.connect(on_save)
    layout.addWidget(btn)
    
    # Keyboard shortcut: Press 'S' to save
    @viewer.bind_key('s')
    def _(viewer):
        on_save()

    # Add control panel to viewer
    viewer.window.add_dock_widget(ctrl, name="Controls", area="right")
    
    # Start napari event loop
    napari.run()

    # ========================================
    # Save Results with DPI Preservation
    # ========================================
    if result['saved']:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        # Extract DPI from original slide
        orig = Image.open(slide_path)
        dpi = orig.info.get('dpi', (72, 72))

        # Save ROI
        Image.fromarray(roi_img).save(out / 'roi.tif', dpi=dpi)
        
        # Save registered slide
        Image.fromarray(result['image']).save(out / 'registered.tif', dpi=dpi)
        
        # Save RGB overlay (R=ROI, G=registered)
        overlay = np.zeros_like(roi_img)
        overlay[:, :, 0] = roi_img[:, :, 0]      # Red channel
        overlay[:, :, 1] = result['image'][:, :, 1]  # Green channel
        Image.fromarray(overlay).save(out / 'overlay.tif', dpi=dpi)

        print(f"✔ 저장 완료: {out}")
        return result['image']

    print("✘ 저장 취소됨.")
    return None

################################################################################
# 3. Processing Functions
################################################################################

def process_files(roi_path, slide_path, output_dir):
    """
    Process a single ROI-slide pair.
    
    Args:
        roi_path: Path to ROI image
        slide_path: Path to slide image
        output_dir: Output directory (subdirectory will be created)
    """
    roi = load_rgb(roi_path)
    slide = load_rgb(slide_path)
    
    if roi is None or slide is None:
        print("❌ 이미지 로드 실패")
        return
    
    # Create subdirectory for this pair
    pair_output = Path(output_dir) / Path(roi_path).stem
    align_and_crop(roi, slide, slide_path, pair_output)


def batch_process(folder_path, output_dir):
    """
    Process all matched ROI-slide pairs in a directory.
    
    Workflow:
        1. Scan folder and match ROI-slide pairs
        2. Process each pair sequentially
        3. User reviews and saves each result interactively
    
    Args:
        folder_path: Directory containing ROI and slide images
        output_dir: Output directory for all results
    """
    pairs = match_pairs(Path(folder_path))
    
    if not pairs:
        print("❌ 매칭되는 파일 쌍이 없습니다.")
        return
    
    for roi, slide in pairs.items():
        print(f"\n>> {roi.name} ↔ {slide.name}")
        process_files(roi, slide, output_dir)
    
    print(f"\n✅ 모두 완료 → {output_dir}")

################################################################################
# 4. Command Line Interface
################################################################################

def main():
    
    parser = argparse.ArgumentParser(
        description="Interactive image registration using napari"
    )
    
    subparsers = parser.add_subparsers(dest='cmd', help='Processing mode')
    
    # Single pair mode
    parser_single = subparsers.add_parser(
        'single', 
        help='Process a single ROI-slide pair'
    )
    parser_single.add_argument(
        '--roi', 
        required=True, 
        help='Path to ROI image (*_XXX_YYY.tif)'
    )
    parser_single.add_argument(
        '--slide', 
        required=True, 
        help='Path to slide image (*.tif)'
    )
    parser_single.add_argument(
        '--out', 
        default='results', 
        help='Output directory (default: results)'
    )
    
    # Batch mode
    parser_batch = subparsers.add_parser(
        'batch', 
        help='Process all pairs in a directory'
    )
    parser_batch.add_argument(
        '--dir', 
        required=True, 
        help='Directory containing ROI and slide images'
    )
    parser_batch.add_argument(
        '--out', 
        default='results', 
        help='Output directory (default: results)'
    )
    
    args = parser.parse_args()
    
    if args.cmd == 'single':
        process_files(args.roi, args.slide, args.out)
    elif args.cmd == 'batch':
        batch_process(args.dir, args.out)
    else:
        parser.print_help()


if __name__ == '__main__':
    sys.exit(main())
