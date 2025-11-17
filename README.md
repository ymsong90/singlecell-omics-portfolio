# Multi-Omics Analysis Portfolio

생물정보학 분석 파이프라인 코드 리뷰 자료입니다. 
질병, 혈액 샘플의 면역 미세환경 분석을 위한 single-cell RNA-seq, imaging mass cytometry, mass cytometry 분석 코드를 포함합니다.

---

## scRNA-seq Analysis Pipeline

mouse PORCN knockout 모델의 위 조직 단일세포 전사체 분석 파이프라인 예시입니다. 
Seurat v5 기반으로 구현했으며, 배치 효과 제거와 세포 타입별 differential gene expression 분석에 중점을 두었습니다.

### 핵심 구현 내용

**데이터 통합 및 배치 보정**
- 4개 샘플(WT, KO 각 2 replicates)의 10X Genomics 데이터 처리
- Harmony를 이용한 배치 효과 제거 (PCA embedding 보정)
- 배치 보정 전후 UMAP 비교를 통한 효과 검증

**세포 타입 어노테이션**
- Canonical marker 기반 주요 세포 타입 분류 (T cell, B cell, myeloid, epithelial, stromal)
- Myeloid 세포 서브셋 재클러스터링 (monocyte/macrophage 서브타입 5개 분리)
- 22개 클러스터 → 11개 주요 세포 타입으로 병합

**Differential expression 분석**
- FindMarkers를 이용한 조건별(WT vs KO) DEG 분석
- 세포 타입별 반복 분석 자동화 (loop 구조)
- Volcano plot 및 heatmap 시각화

**기술적 특징**
- Seurat v5 layer 시스템 처리 (JoinLayers 적용)
- 대용량 데이터 메모리 효율적 처리
- 재현성을 위한 seed 고정 및 파라미터 문서화

---

## Imaging Mass Cytometry Analysis

피부 조직의 spatial proteomics 분석 파이프라인입니다. 
Steinbock 기반 전처리와 Python 기반 공간 분석을 결합하여 세포 간 상호작용 및 조직 구조 분석을 수행했습니다.

### 핵심 구현 내용

**이미지 전처리 및 세포 분할**
- Steinbock를 이용한 IMC raw 데이터 전처리
- mesmer 기반 세포 분할(cell segmentation)
- 40+ 마커 패널의 single-cell 발현 데이터 추출

**공간 분석(Spatial analysis)**
- 세포 이웃 분석(spatial neighborhood analysis)
- 세포 타입 간 거리 계산 및 co-localization 분석
- Melanin 발현 세포의 공간 분포 패턴 정량화 (IMC-melanin staining slide image 정합 분석 예시)

**멜라닌 검출 알고리즘 구현**
- Fontana-masson 염색 슬라이드 이미지 기반 멜라닌 검출 파이프라인 개발
- Napari GUI 기반 interactive annotation 도구 구현
- 멜라닌 양성 세포의 in-depth profiling

**시각화**
- Spatial heatmap (세포 밀도, 마커 발현)
- Neighborhood enrichment plot (세포 간 상호작용)
- Overlay 이미지 (마커 + 세포 경계 + 조직 구조)

**기술적 특징**
- TIF 이미지 메타데이터 처리 (OME-TIFF)
- 대용량 멀티채널 이미지 효율적 로딩 (Zarr 포맷 활용)

---

## Mass Cytometry (CyTOF) Analysis

방광암 환자의 말초혈액 면역세포 프로파일링 파이프라인 예시입니다. 
multiplex barcoding 데이터의 de-barcoding부터 batch correction, cell type annotation, 통계 분석까지 전체 워크플로우를 구현했습니다.

### 핵심 구현 내용

**Multiplex 데이터 처리**
- CATALYST를 이용한 6-plex Live-cell barcode 디바코딩
- DVS bead 기반 신호 정규화 (시간에 따른 signal drift 보정)
- Population-specific cutoff를 통한 정확한 샘플 분리

**배치 통합 및 보정**
- 배치별 Y89 채널 처리 (Batch A: barcode, Batch B: marker)
- 공통 마커 기반 배치 병합
- Harmony를 이용한 PCA embedding 보정

**클러스터링 및 어노테이션**
- FlowSOM self-organizing map 기반 클러스터링 (100 nodes, maxK=50)
- 50개 메타클러스터 → 13개 면역세포 타입 어노테이션
- Monocyte 서브셋 재클러스터링 (5개 서브타입)

**Differential analysis**
- diffcyt를 이용한 differential abundance (DA) 분석
- Differential state (DS) 분석으로 마커 발현 변화 검출
- BCG 치료 반응군 vs 비반응군 비교
- EdgeR (DA) 및 limma (DS) 기반 통계 검정

**시각화**
- t-SNE/UMAP 차원 축소 (Harmony-corrected embeddings)
- Expression heatmap (median marker intensity)
- Stacked barplot with hierarchical clustering
- Condition-specific marker expression comparison

**기술적 특징**
- SingleCellExperiment object 활용
- 배치별 패널 차이 처리 (robust feature matching)
- 메타데이터 기반 자동화된 분석 파이프라인
- 논문 수준의 publication-quality 그림 생성

---

## 기술 스택

**R Packages**
- Single-cell: Seurat (v5), CATALYST, SingleCellExperiment
- Batch correction: harmony, batchelor
- Statistics: diffcyt, edgeR, limma, DESeq2
- Visualization: ggplot2, ComplexHeatmap, patchwork

**Python Libraries**
- Image processing: steinbock, scikit-image
- Spatial analysis: Lisaclust, CATALYST
- Deep learning: deepcell(mesmer)
- Visualization: napari, matplotlib, seaborn

**Workflow Management**
- 모듈화된 스크립트 구조 (01_, 02_, 03_...)
- RData/h5ad 기반 중간 결과 저장
- 파라미터 문서화 및 재현 가능한 분석

---

## 분석 특징

이 포트폴리오의 코드는 실제 연구 프로젝트에서 사용된 분석 파이프라인을 기반으로 합니다. 각 분석은 다음과 같은 공통적인 특징을 가집니다:

- **재현성**: 모든 random seed 고정, 파라미터 명시적 문서화
- **효율성**: 대용량 데이터 처리를 위한 메모리 최적화
- **확장성**: 모듈화된 구조로 새로운 데이터셋 적용 용이

---

## Contact

**송영민 (Youngmin Song)**  
M.S. Candidate in Clinical Medical Sciences  
College of Medicine, Seoul National University  
Seoul National University Hospital
