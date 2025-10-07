# NGS中心型91病原体検査システム 最適化戦略 v2.0
## PCR依存を最小化したメタゲノムNGS主導アプローチ

**作成日**: 2025年10月8日
**改訂**: v2.0 - PCR minimal strategy
**対象**: 異種移植用ドナーブタ（Yucatan miniature pig, 3KO-7TG-59PERV）
**サンプル**: 血漿由来 cfDNA/cfRNA

---

## エグゼクティブサマリー

### 戦略転換: **NGS Single Platform Strategy**

```
【旧戦略 v1.0】
NGS (主) + マルチプレックスPCR (確認) + 血清学 (補助)
└─ 問題点:
   ├─ PCR プライマー設計コスト: 2-3万円 × 91病原体
   ├─ 変異対応の継続的コスト: 年間50-100万円
   ├─ 新規変異株への対応遅延: 2-4週間
   └─ 偽陰性リスク（プライマーミスマッチ）

【新戦略 v2.0】
Deep NGS Metagenomics (100%) + Computational Validation
└─ 利点:
   ✓ プライマー設計不要 → コスト削減
   ✓ 変異株・亜種を自動検出
   ✓ リアルタイム系統解析
   ✓ 単一プラットフォームで完結
   ✓ 定量精度向上（リード数ベース）
```

### コスト比較

| 項目 | v1.0 (PCR併用) | v2.0 (NGS単独) | 削減額 |
|------|----------------|----------------|--------|
| **Phase 1** | 150,000円 + PCR 30,000円 | **120,000円** | **-60,000円** |
| **Phase 2** | 50,000円 + PCR 20,000円 | **35,000円** | **-35,000円** |
| **Phase 3** | 35,000円 + PCR 15,000円 | **25,000円** | **-25,000円** |
| **継続コスト** | プライマー更新 年50万円 | **データベース更新 年10万円** | **-40万円/年** |

**結論**: NGS単独戦略により**30-60%のコスト削減** + **変異対応の自動化**

---

## 1. PCR排除の技術的根拠

### 1.1 PCRの問題点（詳細分析）

```
【Problem 1: プライマー設計の複雑性とコスト】

初期設計コスト:
├─ 91病原体 × 各3-5領域 = 273-455 プライマーセット
├─ 設計・合成: 10,000-15,000円/セット
├─ バリデーション: 20,000円/セット
└─ 総額: 8,190,000-13,650,000円

継続的メンテナンスコスト（年間）:
├─ 変異モニタリング: 200,000円
├─ プライマー再設計: 15セット × 30,000円 = 450,000円
├─ クロスバリデーション: 200,000円
└─ 年間合計: 850,000円

【Problem 2: 変異による偽陰性リスク】

プライマー結合部位の変異頻度:
├─ RNAウイルス（インフルエンザ等）: 10^-3 - 10^-4 /塩基/複製
├─ レトロウイルス（PERV）: 10^-5 /塩基/複製
├─ DNAウイルス: 10^-6 - 10^-8 /塩基/複製
├─ 細菌: 10^-9 - 10^-10 /塩基/複製

プライマー（20-25塩基）での変異確率:
├─ RNAウイルス: 2-5%/世代 → 高頻度の偽陰性
├─ PERV: 0.02-0.05%/世代 → 中程度のリスク
└─ 細菌: 極めて低い

実例:
├─ ブタインフルエンザ: 亜型間でプライマー非互換性
├─ PRRSV: 北米型・欧州型で異なるプライマー必要
├─ PERV: A/B/C各型 + リコンビナントで複雑化
└─ サルモネラ: 2,500以上の血清型

【Problem 3: マルチプレックスの技術的制約】

同時検出可能数:
├─ TaqMan方式: 最大4-5色 → 5病原体/ウェル
├─ 91病原体 ÷ 5 = 18-20反応必要
├─ クロストーク・干渉のリスク
└─ 最適化の複雑性とコスト

【Problem 4: 定量精度の限界】

PCRの定量問題:
├─ 増幅効率の変動: 85-105%/サイクル
├─ プラトー効果: 高濃度サンプルで非線形
├─ 阻害物質の影響: 血漿中のヘム、免疫グロブリン等
└─ スタンダードカーブの維持コスト

対してNGSの定量:
✓ リード数 = 直接的な分子数反映
✓ 線形性: 6-7 log dynamic range
✓ 複数病原体の同時定量
✓ スタンダード不要（spike-in使用）
```

### 1.2 NGS Deep Sequencingの優位性

```
【Advantage 1: 変異株への完全対応】

プライマーフリー検出:
├─ ライブラリ調製はランダムフラグメント化
├─ アダプター配列のみで増幅（病原体配列非依存）
├─ 配列全体を取得 → 変異部位も含めて検出
└─ 新規変異株を自動的に検出

系統解析の自動化:
├─ 全ゲノム/トランスクリプトーム配列取得
├─ リアルタイムで系統樹作成
├─ 既知株との相同性解析
└─ 薬剤耐性変異の検出

【Advantage 2: 定量精度と再現性】

デジタル定量:
├─ 各病原体のリード数 = 分子数
├─ 正規化: TPM (Transcripts Per Million) or RPM
├─ UMI (Unique Molecular Identifier) 使用で絶対定量
└─ CV (変動係数) < 15% （PCRの20-30%と比較）

ダイナミックレンジ:
├─ 10^1 - 10^7 copies/mL を同時検出
├─ 希少病原体と優勢病原体を同一ラン検出
└─ リード深度調整で感度制御

【Advantage 3: 多検体並列処理】

バーコーディング:
├─ 最大384検体/ラン（NovaSeq 6000）
├─ 96検体/ラン（NextSeq 2000）
├─ 24検体/ラン（MiSeq）
└─ 検体数増加でコスト逓減

【Advantage 4: データの二次利用価値】

蓄積データの活用:
├─ 長期疫学モニタリング
├─ 新規病原体の遡及検出
├─ AI/機械学習モデルのトレーニング
├─ メタゲノムデータベース構築
└─ 論文発表・知的財産化
```

---

## 2. NGS単独戦略の詳細設計

### 2.1 シーケンス深度最適化

#### 2.1.1 必要リード数の理論計算

```
【前提条件】
- サンプル: 5mL 血漿
- cfDNA/RNA total: 50-200ng (平均100ng)
- ホスト除去効率: 95%
- 病原体核酸比率: 5% (残り)
- 検出したい最低コピー数: 10 copies/mL (LOD)

【計算】
血漿中総核酸分子数:
├─ 100ng ≈ 3×10^10 分子 (仮定: 平均150bp fragment)
├─ 病原体分子数 (5%): 1.5×10^9 分子
└─ 最低検出コピー (10 copies/mL): 50 copies in 5mL

必要リード深度:
├─ 検出確率 >99% (ポアソン分布)
├─ λ = 期待リード数 ≥ 10 reads
├─ 病原体比率 = 50 copies / 1.5×10^9 total = 3.3×10^-8
├─ 必要総リード数 = 10 / 3.3×10^-8 = 3×10^8 reads
└─ 安全マージン ×2 = 6×10^8 reads

【実用的設定】
Phase 1 (外注): 50-100M reads/sample
└─ 余裕を持った検出感度
└─ 未知病原体検出にも対応

Phase 2-3 (内製): 30-50M reads/sample
└─ コスト最適化
└─ 既知病原体フォーカス

Phase 4 (最適化): 20-30M reads/sample
└─ UMI使用で効率化
└─ ターゲットエンリッチメント併用
```

#### 2.1.2 リード深度別コストと感度

| リード数 | 検出LOD | コスト/検体 | 推奨フェーズ | 備考 |
|---------|---------|------------|-------------|------|
| **10M reads** | 100 copies/mL | 60,000円 | - | 感度不足 |
| **20M reads** | 50 copies/mL | 80,000円 | Phase 4 | 最適化後 |
| **30M reads** | 20 copies/mL | 90,000円 | Phase 3 | 標準 |
| **50M reads** | 10 copies/mL | 120,000円 | **Phase 1-2** | **推奨** |
| **100M reads** | 5 copies/mL | 180,000円 | 研究用 | 過剰 |

**結論**: **Phase 1-2は50M reads、Phase 3以降は30M readsが最適**

### 2.2 ホスト遺伝子除去の最適化

#### 2.2.1 ホスト除去の重要性

```
【ホスト除去なしの場合】
血漿cfDNA/RNA組成:
├─ ブタ宿主由来: 99%
├─ 病原体由来: 0.01-1%
└─ その他（食餌、環境）: 0.01-0.1%

50M readsの配分:
├─ ホスト: 49.5M reads (99%) → 無駄
├─ 病原体: 50,000-500,000 reads (0.1-1%)
└─ ノイズ: 50,000 reads

【ホスト除去 95%後】
50M readsの配分:
├─ ホスト残存: 2.5M reads (5%)
├─ 病原体: 950,000-9,500,000 reads (19-95%) → 有効リード
└─ 19-190倍の感度向上！
```

#### 2.2.2 多段階ホスト除去戦略

```yaml
STRATEGY: "Triple-Layer Host Depletion"

Layer_1_Enzymatic_Depletion:
  Method: "CpG-methylation based"
  Kit: "NEBNext Microbiome DNA Enrichment Kit"
  Target: "Methylated host DNA (CpG islands)"
  Efficiency: "70-85%"
  Cost: "5,000円/sample"
  Time: "2 hours"

Layer_2_Hybridization_Capture:
  Method: "Biotin-streptavidin depletion"
  Kit: "QIAseq FastSelect (Custom pig probes)"
  Target:
    - "Sus scrofa ribosomal RNA (99% of RNA)"
    - "Highly expressed pig genes (top 1000)"
    - "Modified pig transgenes (3KO-7TG)"
  Efficiency: "90-98%"
  Cost: "8,000円/sample"
  Time: "4 hours"
  Notes: "RNA特異的、DNA除去には併用不可"

Layer_3_Computational_Depletion:
  Method: "In-silico host read removal"
  Tool: "BWA-MEM2 + Bowtie2 (ultra-sensitive mode)"
  Reference:
    - "Sus_scrofa_11.1 genome"
    - "Custom modified pig genome (3KO-7TG-59PERV)"
    - "Pig repeat sequences (RepeatMasker)"
  Efficiency: "95-99.5% (cumulative)"
  Cost: "計算コストのみ（1,000円/sample）"
  Time: "1-2 hours"

COMBINED_EFFICIENCY:
  DNA_pathway:
    - Layer 1 (enzymatic): 80% removal
    - Layer 3 (computational): 98% of remaining
    - Total: 99.6% removal
    - Pathogen enrichment: 250x

  RNA_pathway:
    - Layer 2 (rRNA depletion): 95% removal
    - Layer 3 (computational): 98% of remaining
    - Total: 99.9% removal
    - Pathogen enrichment: 1000x

COST_BENEFIT:
  Total_cost: "14,000円/sample (Layer 1+2)"
  Benefit: "200-1000倍の感度向上"
  ROI: "極めて高い（必須投資）"
```

### 2.3 ライブラリ調製の最適化

#### 2.3.1 Dual Library Strategy

```yaml
LIBRARY_PREP_WORKFLOW:

DNA_Library:
  Input: "cfDNA (post-host-depletion, 5-20ng)"
  Method: "Tagmentation-based"
  Kit: "Illumina DNA Prep (旧 Nextera Flex)"

  Advantages:
    - "Low input requirement (1ng可能)"
    - "Fast workflow (90 min)"
    - "Uniform coverage"
    - "UMI compatible"

  Steps:
    1_Tagmentation:
      - Transposase cuts DNA + adds adapters
      - Fragment size: 150-300bp peak
      - Time: 15 min

    2_PCR_Amplification:
      - Cycles: 5-7 (minimal amplification bias)
      - Dual indexing for multiplexing
      - Time: 45 min

    3_Cleanup:
      - SPRI beads (0.8x ratio)
      - Time: 30 min

  QC_Metrics:
    - Concentration: >2 nM
    - Size: 200-400bp
    - GC bias: <10%

  Cost: "8,000円/sample"

RNA_Library:
  Input: "cfRNA (post-rRNA-depletion, 10-50ng)"
  Method: "Random priming + strand-specific"
  Kit: "Illumina Stranded Total RNA Prep"

  Advantages:
    - "Detects all RNA (mRNA + viral RNA)"
    - "Strand information preserved"
    - "Quantitative (no 3' bias)"

  Steps:
    1_rRNA_Depletion:
      - QIAseq FastSelect (pig-specific probes)
      - Removes >95% rRNA
      - Time: 4 hours

    2_Fragmentation:
      - Chemical fragmentation (Mg2+, heat)
      - Target: 150-200bp
      - Time: 15 min

    3_Reverse_Transcription:
      - Random hexamer priming
      - dUTP 2nd strand synthesis (strand-specific)
      - Time: 1 hour

    4_PCR_Amplification:
      - Cycles: 12-15
      - Dual indexing
      - Time: 45 min

    5_Cleanup:
      - SPRI beads (0.9x ratio)
      - Time: 30 min

  QC_Metrics:
    - Concentration: >2 nM
    - Size: 250-450bp
    - rRNA remaining: <5%
    - Strand specificity: >90%

  Cost: "12,000円/sample"

POOLING_STRATEGY:
  Ratio: "DNA:RNA = 1:1 (molar ratio)"
  Reason: "Balanced detection of DNA/RNA pathogens"

  Alternative_for_cost_reduction:
    DNA_only:
      - "DNAライブラリのみ（20,000円削減）"
      - "RNAウイルスは逆転写後に検出"
      - "感度低下: 10-50%"
      - "推奨: Noいーーーーーーー（RNAウイルス重要）"

TOTAL_LIBRARY_COST: "20,000円/sample"
```

#### 2.3.2 UMI (Unique Molecular Identifier) の導入

```yaml
UMI_IMPLEMENTATION:
  Timing: "Phase 3以降で導入"

  Purpose:
    - "PCR duplication除去"
    - "絶対定量の実現"
    - "低頻度変異の検出"

  Method:
    - "12bp random barcode at adapter ligation"
    - "Pre-PCR molecular tagging"
    - "Post-sequencing deduplication"

  Benefits:
    Quantification_accuracy:
      - Without_UMI: "CV = 20-30%"
      - With_UMI: "CV = 5-10%"
      - Improvement: "2-6倍の精度向上"

    Sensitivity:
      - "低コピー数病原体の検出感度向上"
      - "LOD改善: 10 copies → 5 copies"

    Variant_detection:
      - "PCR errorとreal variantの区別"
      - "0.1%の変異検出可能"

  Cost:
    Additional: "+3,000円/sample (UMI adapter)"
    Computational: "+500円/sample (deduplication)"
    Total: "+3,500円/sample"

  ROI_Analysis:
    Break_even: "Phase 3以降（検体数増加時）"
    Recommendation: "Phase 3で標準実装"
```

### 2.4 シーケンシングプラットフォーム選択

#### 2.4.1 Phase別最適プラットフォーム

```yaml
PHASE_1_External_Sequencing:
  Platform: "NovaSeq 6000 (外注)"
  Provider:
    - "タカラバイオ"
    - "Macrogen"
    - "GENEWIZ"

  Specifications:
    Flowcell: "SP (Single lane)"
    Output: "650-800M reads/lane"
    Multiplexing: "12-16 samples/lane"
    Reads_per_sample: "40-66M reads"
    Read_length: "PE150"

  Cost_breakdown:
    Sequencing: "400,000円/lane"
    Per_sample: "25,000-33,000円"
    Library_prep: "20,000円/sample"
    Total: "45,000-53,000円/sample"

  Turnaround: "2-3 weeks"

  Advantages:
    ✓ "初期投資不要"
    ✓ "最新装置を利用可能"
    ✓ "技術習得の猶予"

  Disadvantages:
    ✗ "ターンアラウンド長い"
    ✗ "データ機密性懸念"
    ✗ "長期的コスト高"

PHASE_2_Partial_Internal:
  Platform: "NextSeq 2000 (共同利用 or リース)"
  Location:
    Option_A: "大学共同利用施設"
    Option_B: "受託会社の時間貸し"
    Option_C: "リース契約"

  Specifications:
    Flowcell: "P2 (200 cycles)"
    Output: "400-1,200M reads"
    Multiplexing: "8-24 samples"
    Reads_per_sample: "16-150M reads"
    Read_length: "PE100 or PE150"

  Cost_breakdown:
    Flowcell: "150,000円/run"
    Per_sample (24 samples): "6,250円"
    Library_prep: "20,000円/sample (内製)"
    Total: "26,250円/sample"

  共同利用の場合:
    Usage_fee: "10,000円/sample + 消耗品"
    Total: "30,000円/sample"

  Turnaround: "1 week"

  Advantages:
    ✓ "コスト削減（Phase 1比 40%削減）"
    ✓ "ターンアラウンド短縮"
    ✓ "NGS技術の習得"

  Disadvantages:
    ✗ "共同利用の予約制約"
    ✗ "完全内製化できない"

PHASE_3_Full_Internal:
  Platform: "MiSeq (購入) or NextSeq 2000 (購入)"

  Decision_criteria:
    MiSeq:
      - Sample_throughput: "<12 samples/month"
      - Budget: "Limited (15M円)"
      - Flexibility: "High (即座に利用可能)"

    NextSeq_2000:
      - Sample_throughput: ">24 samples/month"
      - Budget: "Sufficient (25M円)"
      - Cost_efficiency: "High (長期的に有利)"

  MiSeq_Option:
    Equipment_cost: "15,000,000円"
    Reagent_cost: "80,000円/run (v3, 600 cycles)"
    Samples_per_run: "24 samples"
    Per_sample_reagent: "3,333円"

    Output:
      - "25M reads/run"
      - "~1M reads/sample (24 samples)"
      - "十分だが余裕少ない"

    Total_per_sample:
      Library: "20,000円"
      Sequencing: "3,333円"
      Computational: "1,000円"
      Total: "24,333円/sample"

    Annual_cost (24 samples):
      Reagents: "80,000円 × 12 runs = 960,000円"
      Maintenance: "1,500,000円/year"
      Personnel: "8,000,000円"
      Total: "10,460,000円"
      Per_sample: "435,833円"

  NextSeq_2000_Option:
    Equipment_cost: "25,000,000円"
    Reagent_cost: "150,000円/run (P2, 200 cycles)"
    Samples_per_run: "24-48 samples"
    Per_sample_reagent: "3,125-6,250円"

    Output:
      - "400-1,200M reads/run"
      - "16-50M reads/sample"
      - "十分な余裕"

    Total_per_sample (24 samples):
      Library: "20,000円"
      Sequencing: "6,250円"
      Computational: "1,500円"
      Total: "27,750円/sample"

    Annual_cost (24 samples):
      Reagents: "150,000円 × 12 runs = 1,800,000円"
      Maintenance: "2,500,000円/year"
      Personnel: "10,000,000円"
      Total: "14,300,000円"
      Per_sample: "595,833円"

    But_scaling (48 samples):
      Per_sample: "327,083円"

  Recommendation:
    Low_throughput: "MiSeq (24検体/年以下)"
    High_throughput: "NextSeq 2000 (24検体/年以上)"
    Future_proof: "NextSeq 2000 (拡張性重視)"

PHASE_4_Optimization:
  Platform: "NextSeq 2000 (最適化運用)"

  Optimization_strategies:
    1_Batch_processing:
      - "48-96 samples/run"
      - "Per-sample cost: 3,125円 (reagent only)"

    2_Targeted_enrichment:
      - "PMDA 91 pathogen panel probes"
      - "Reduces required depth: 30M → 10M"
      - "Cost reduction: 40%"

    3_AI_predictive_QC:
      - "リアルタイム品質モニタリング"
      - "Failed run early detection"
      - "Reagent waste reduction: 10%"

  Target_cost: "20,000円/sample (all-in)"
```

#### 2.4.2 コスト比較サマリー

| Phase | Platform | Reagent | 総コスト/検体 | TAT |
|-------|----------|---------|--------------|-----|
| **1** | NovaSeq (外注) | 含まれる | **50,000円** | 2-3週 |
| **2** | NextSeq (共同利用) | 6,250円 | **27,000円** | 1週 |
| **3a** | MiSeq (内製) | 3,333円 | **24,000円** | 3-4日 |
| **3b** | NextSeq (内製24検体) | 6,250円 | **28,000円** | 3-4日 |
| **3c** | NextSeq (内製48検体) | 3,125円 | **24,500円** | 3-4日 |
| **4** | NextSeq (最適化) | 2,000円 | **20,000円** | 2-3日 |

---

## 3. バイオインフォマティクスパイプライン（PCRフリー版）

### 3.1 パイプライン全体設計

```yaml
PIPELINE: "PMDA_Pathogen_NGS_Pipeline_v2.0"
PHILOSOPHY: "PCR-free, NGS-centric, Computationally intensive"

ARCHITECTURE:
  Input: "Paired-end FASTQ files (DNA + RNA)"
  Output: "Pathogen detection report + Quantification + Variants"
  Language: "Nextflow DSL2 (workflow) + Python (scripts)"
  Containerization: "Docker/Singularity (reproducibility)"

WORKFLOW_STAGES:

  #================================
  # STAGE 1: Quality Control
  #================================
  Stage_1_QC:
    Tools:
      - FastQC: "v0.12.1"
      - MultiQC: "v1.14"
      - fastp: "v0.23.4 (trimming + filtering)"

    Steps:
      1_Raw_QC:
        - "Raw read quality assessment"
        - "Adapter detection"
        - "GC content analysis"
        - "Duplication rate estimation"

      2_Trimming:
        Tool: "fastp"
        Parameters:
          - "--cut_front --cut_tail"  # Quality-based trimming
          - "--cut_mean_quality 20"   # Q20 threshold
          - "--length_required 50"    # Minimum length
          - "--detect_adapter_for_pe" # Auto-detect adapters
          - "--dedup"                 # Remove exact duplicates (if no UMI)

      3_Post_trim_QC:
        - "Verify Q30 >85%"
        - "Check read length distribution"
        - "Assess data loss <10%"

    Quality_Gates:
      PASS_criteria:
        - Q30_score: ">85%"
        - Total_reads: ">20M"
        - Mean_length: ">100bp"
        - Duplication: "<40%"

      FAIL_action:
        - "Flag for re-sequencing"
        - "Notify operator"

    Output: "Trimmed, filtered FASTQ"
    Time: "15-30 min"
    CPU: "8 cores"
    RAM: "16GB"

  #================================
  # STAGE 2: Host Removal
  #================================
  Stage_2_Host_Removal:
    Reference_Genomes:
      Sus_scrofa:
        Version: "Sscrofa11.1"
        Source: "Ensembl"
        Size: "2.5GB"
        Chromosomes: "18 + X + Y + MT"

      Modified_Pig:
        Description: "3KO-7TG-59PERV custom genome"
        Base: "Sscrofa11.1"
        Modifications:
          Knockouts:
            - "GGTA1 (chr6)"
            - "CMAH (chr7)"
            - "B4GALNT2 (chr2)"
          Transgenes:
            - "hCD46, hCD55, hTHBD, hPROCR"
            - "hCD47, hTNFAIP3, hHMOX1"
          PERV_inactivations: "59 loci (distributed)"

        Build_method:
          - "Custom FASTA with edits"
          - "CRISPR edit positions annotated"
          - "Transgene sequences inserted"

      Combined_index:
        Tool: "BWA-MEM2 index + Bowtie2 index"
        Size: "8GB (indexed)"

    Alignment_Strategy:
      Primary_alignment:
        Tool: "BWA-MEM2"
        Mode: "Very-sensitive"
        Parameters:
          - "-k 19"          # Seed length
          - "-w 100"        # Band width
          - "-d 100"        # Z-dropoff
          - "-r 1.5"        # Re-seed trigger
          - "-A 1 -B 4"     # Match/mismatch scores
          - "-O 6 -E 1"     # Gap penalties
          - "-T 30"         # Minimum score (stringent)

        Purpose: "Remove perfect and near-perfect host matches"
        Expected_removal: "85-95%"

      Secondary_alignment:
        Tool: "Bowtie2"
        Mode: "Sensitive-local"
        Parameters:
          - "--very-sensitive-local"
          - "--score-min L,-0.6,-0.6"
          - "--no-mixed --no-discordant"

        Purpose: "Remove partial host matches (fragmented reads)"
        Expected_additional_removal: "3-5%"

      Combined_efficiency: "95-99%"

    Extract_non_host:
      Tool: "samtools view -f 4"  # Unmapped reads
      Format: "FASTQ (non-host enriched)"
      Expected_output: "1-5% of input reads"

    Quality_check:
      - "Host removal rate >95%"
      - "Pathogen spike-in recovery >90%"
      - "Non-specific loss <5%"

    Output: "Non-host FASTQ (pathogen-enriched)"
    Time: "1-2 hours"
    CPU: "32 cores"
    RAM: "64GB"
    Storage: "Temporary SAM/BAM: 50-100GB"

  #================================
  # STAGE 3: Pathogen Detection (Multi-module)
  #================================
  Stage_3_Pathogen_Detection:

    #----- Module A: Rapid Taxonomic Classification -----
    Module_A_Kraken2_Bracken:
      Tool: "Kraken2 v2.1.2 + Bracken v2.8"

      Database:
        Name: "PMDA_Custom_Kraken2_DB"
        Components:
          - "Kraken2 Standard DB (archaea, bacteria, viruses, plasmids)"
          - "PMDA 91 designated pathogens (augmented)"
          - "Viral RefSeq (complete genomes)"
          - "Custom PERV sequences (A, B, C, recombinants)"

        Size: "70GB (indexed)"
        Update: "Monthly (automated script)"

      Workflow:
        1_Kraken2_classification:
          Command: >
            kraken2 --db PMDA_DB
                    --paired read1.fq read2.fq
                    --threads 32
                    --confidence 0.1
                    --minimum-base-quality 20
                    --report kraken2_report.txt
                    --output kraken2_output.txt

          Output:
            - "Read-level classification"
            - "Taxonomic abundance report"

          Time: "5-10 min"

        2_Bracken_quantification:
          Command: >
            bracken -d PMDA_DB
                    -i kraken2_report.txt
                    -o bracken_output.txt
                    -r 150
                    -l S  # Species level
                    -t 10 # Minimum reads

          Output:
            - "Re-estimated abundance (corrected for genome size)"
            - "Species-level quantification"

          Time: "1-2 min"

      Sensitivity:
        - "Detects pathogens with >10 reads"
        - "Genus/Species level classification"
        - "Low computational cost"

      Limitations:
        - "Database-dependent (known pathogens only)"
        - "May misclassify closely related species"
        - "Not suitable for novel pathogen discovery"

      Use_case: "First-pass screening for known PMDA pathogens"

    #----- Module B: Virus-Specific Detection -----
    Module_B_RVDB_BLAST:
      Tool: "DIAMOND v2.1 (BLASTX mode, 1000x faster than BLAST)"

      Database:
        Name: "RVDB v30.0 (Reference Viral Database)"
        Source: "https://rvdb.dbi.udel.edu/"
        Description: "Manually curated viral sequences"

        Contents:
          - "Virus protein sequences: 3.5M"
          - "Virus nucleotide sequences: 2.8M"
          - "PMDA viral pathogens: 35 species (all strains)"

        Customization:
          - "Added: PERV variants (A, B, C, A/C recombinant)"
          - "Added: Pig-specific viruses (non-PMDA)"
          - "Removed: Plant/marine viruses (reduce noise)"

        Size: "25GB (DIAMOND indexed)"
        Update: "Biannually"

      Workflow:
        1_DIAMOND_BLASTX:
          Command: >
            diamond blastx
                    --query non_host.fasta
                    --db RVDB_diamond.dmnd
                    --threads 32
                    --evalue 1e-5
                    --max-target-seqs 5
                    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids
                    --out diamond_output.txt

          Parameters:
            - E-value: "1e-5 (stringent)"
            - Identity: "Auto (no hard cutoff)"
            - Max_targets: "5 (top 5 hits per read)"

          Time: "30-60 min"

        2_LCA_assignment:
          Tool: "Custom Python script"
          Method: "Lowest Common Ancestor (LCA) algorithm"

          Logic:
            - "If top hits agree → assign to species"
            - "If top hits conflict → assign to genus/family"
            - "Filter ambiguous assignments (>5 taxa)"

          Output: "Virus taxonomy table with read counts"
          Time: "5 min"

      Sensitivity:
        - "Detects divergent viral strains (60-70% identity OK)"
        - "Protein-level homology (more sensitive than nucleotide)"
        - "Cross-species detection"

      Specificity:
        - "High (E-value 1e-5)"
        - "LCA reduces false positives"

      Use_case: "Comprehensive viral screening including variants"

    #----- Module C: Bacterial/Fungal Detection -----
    Module_C_MetaPhlAn4:
      Tool: "MetaPhlAn v4.0"

      Database:
        Name: "mpa_vJan21_CHOCOPhlAnSGB_202103"
        Contents:
          - "1.1M unique clade-specific markers"
          - "Bacteria, Archaea, Eukaryotes, Viruses"
          - "Species Genome Bins (SGBs): 26,970"

        Augmentation:
          - "Added: PMDA bacterial pathogens (27 species)"
          - "Added: PMDA fungal pathogens (2 species)"
          - "Custom markers for Salmonella, Brucella, Leptospira"

        Size: "17GB"

      Workflow:
        Command: >
          metaphlan non_host.fastq.gz
                   --bowtie2db mpa_vJan21
                   --nproc 32
                   --input_type fastq
                   --tax_lev s  # Species level
                   --min_cu_len 2000
                   --stat_q 0.2
                   --output metaphlan_profile.txt
                   --CAMI_format_output

        Output:
          - "Taxonomic profiling (relative abundance %)"
          - "Species-level resolution"
          - "CAMI format (standardized)"

        Time: "20-30 min"

      Sensitivity:
        - "Marker-based approach (high specificity)"
        - "Quantitative (relative abundance)"
        - "Strain-level resolution (some species)"

      Use_case: "Bacterial/Fungal quantification and profiling"

    #----- Module D: Parasite Detection -----
    Module_D_Parasite_BLAST:
      Tool: "BLASTN v2.13 (nucleotide-nucleotide)"

      Database:
        Name: "PMDA_Parasite_DB"
        Sources:
          - "EuPathDB (Eukaryotic Pathogen Database)"
          - "NCBI Parasite genomes"
          - "PMDA designated parasites: 19 species"

        Contents:
          Protozoa:
            - "Toxoplasma gondii"
            - "Cryptosporidium spp."
            - "Balantidium coli"
            - "Babesia spp."
            - "Trypanosoma spp."
            - "Coccidia, Sarcocystis"

          Helminths:
            - "Ascaris suum (pig roundworm)"
            - "Toxocara spp."
            - "Echinococcus spp."
            - "Trichuris suis (whipworm)"
            - "Strongyloides spp."
            - "Taenia solium (pork tapeworm)"
            - "Trichinella spiralis"
            - "Other (15 species total)"

        Size: "45GB (BLAST formatted)"

      Workflow:
        Command: >
          blastn -query non_host.fasta
                 -db PMDA_Parasite_DB
                 -num_threads 32
                 -evalue 1e-10
                 -max_target_seqs 5
                 -outfmt '6 std staxids'
                 -out blastn_parasite.txt

        Filtering:
          - "E-value < 1e-10 (very stringent)"
          - "Identity >90% (high confidence)"
          - "Alignment length >100bp"

        Time: "1-2 hours"

      Challenge:
        - "Parasites have large genomes (100M-1Gb)"
        - "cfDNA/RNA detection may be low sensitivity"
        - "Environmental contamination risk"

      Validation:
        - "Cross-reference with microscopy (if available)"
        - "Require >100 reads for positive call"

      Use_case: "Parasite screening (confirmatory method recommended)"

    #----- Module E: PERV Monitoring (Critical) -----
    Module_E_PERV_Analysis:
      Importance: "最重要 - PMDA特別監視項目"

      Strategy: "Multi-tiered PERV detection + quantification"

      Tier_1_Genome_mapping:
        Tool: "BWA-MEM2 + minimap2"
        Reference:
          - "PERV-A consensus (GenBank: AF435967)"
          - "PERV-B consensus (GenBank: AF435968)"
          - "PERV-C consensus (GenBank: AF435969)"
          - "PERV-A/C recombinant (env region)"

        Mapping:
          Command: >
            bwa mem -t 32
                    -k 19 -w 100
                    -M  # Mark shorter splits as secondary
                    PERV_ref.fa
                    non_host_R1.fq non_host_R2.fq
                    | samtools view -bS -F 2308 -
                    | samtools sort -o PERV_mapped.bam

          Filter: "Primary alignments only, MAPQ >20"

        Output:
          - "PERV-A mapped reads"
          - "PERV-B mapped reads"
          - "PERV-C mapped reads"
          - "Recombinant candidates"

        Time: "15 min"

      Tier_2_Variant_calling:
        Tool: "FreeBayes v1.3 (bayesian variant caller)"

        Purpose:
          - "Detect strain variations"
          - "Identify recombination breakpoints"
          - "Monitor evolution"

        Command: >
          freebayes -f PERV_ref.fa
                    -p 1  # Haploid (viral)
                    --min-alternate-fraction 0.01  # 1% sensitivity
                    --min-coverage 10
                    PERV_mapped.bam
                    > PERV_variants.vcf

        Analysis:
          - "Envelope (env) gene mutations"
          - "Receptor binding domain changes"
          - "Drug resistance markers (if applicable)"

        Time: "10 min"

      Tier_3_Quantification:
        Method: "Read count normalization"

        Formula:
          PERV_copies_per_mL = (PERV_reads / Total_reads) × (Total_input_molecules / Sample_volume_mL)

        Normalization:
          - "TPM (Transcripts Per Million)"
          - "RPKM (Reads Per Kilobase per Million)"
          - "Spike-in control (ERCC RNA)"

        Thresholds:
          PERV_A:
            - Detection: ">10 reads"
            - Quantification: ">100 reads"
            - Alert: ">1,000 copies/mL"

          PERV_B:
            - Detection: ">10 reads"
            - Quantification: ">100 reads"
            - Alert: ">1,000 copies/mL"

          PERV_C:
            - Detection: ">5 reads (rare)"
            - Alert: "Any detection (要確認)"

          Recombinant:
            - Detection: "Breakpoint evidence"
            - Alert: "Immediate notification"

      Tier_4_Phylogenetic_analysis:
        Tool: "IQ-TREE v2.2 (maximum likelihood)"

        Workflow:
          1_Consensus_generation:
            - "Generate consensus from mapped reads"
            - "Minimum coverage: 10x"

          2_Multiple_alignment:
            - "MAFFT v7.505 (--auto mode)"
            - "Align to known PERV references"

          3_Tree_building:
            - "IQ-TREE with ModelFinder"
            - "1000 bootstrap replicates"

          4_Interpretation:
            - "Identify closest reference strain"
            - "Detect recombination events"
            - "Track temporal evolution (longitudinal samples)"

        Time: "30 min"

      Reporting:
        - "PERV type (A/B/C/recombinant)"
        - "Viral load (copies/mL)"
        - "Key mutations"
        - "Phylogenetic placement"
        - "Risk assessment (低/中/高)"

      Alert_system:
        High_priority:
          - "PERV-C detection"
          - "Novel recombinant"
          - "Viral load >10,000 copies/mL"
          - "High-risk mutations"

        Action:
          - "Immediate email notification"
          - "Automatic report generation"
          - "Recommend confirmatory testing"

    #----- Module F: Novel Pathogen Discovery -----
    Module_F_De_Novo_Assembly:
      Purpose: "未知・新興病原体の検出"

      Strategy: "Assembly-based detection"

      Workflow:
        1_De_novo_assembly:
          Tool: "metaSPAdes v3.15 (metagenomic mode)"

          Command: >
            metaspades.py
                    -1 non_host_R1.fq
                    -2 non_host_R2.fq
                    -k 21,33,55,77
                    --threads 32
                    --memory 250
                    -o metaspades_output

          Parameters:
            - K-mer sizes: "21, 33, 55, 77 (multi-k strategy)"
            - Memory: "250GB"
            - Expected_contigs: "10,000-100,000"

          Output: "Assembled contigs (FASTA)"
          Time: "4-8 hours"

        2_Contig_filtering:
          Criteria:
            - Length: ">500bp"
            - Coverage: ">5x"
            - Non-host: "Re-screen against pig genome"

          Expected_filtered: "1,000-5,000 contigs"

        3_Taxonomic_annotation:
          Tool: "DIAMOND BLASTX against nr database"

          Command: >
            diamond blastx
                    --query filtered_contigs.fa
                    --db nr.dmnd
                    --threads 32
                    --evalue 1e-5
                    --max-target-seqs 10
                    --outfmt 6
                    --out contigs_nr_blast.txt

          Time: "2-4 hours (nr database is huge: 300GB)"

        4_Novel_pathogen_detection:
          Criteria:
            - "Viral contigs: No match in RVDB"
            - "Bacterial contigs: Novel species (16S rRNA divergence >3%)"
            - "Eukaryotic contigs: No match in known parasites"

          Validation:
            - "Manual review of top candidates"
            - "BLAST against nt/nr (comprehensive)"
            - "Coverage plot analysis (even coverage = real)"
            - "Gene prediction (ORF finding)"

          Reporting:
            - "Novel pathogen candidates (ranked by evidence)"
            - "Taxonomic assignment (best guess)"
            - "Recommendation for follow-up (PCR, Sanger)"

      Limitations:
        - "Computationally expensive (8+ hours)"
        - "High false positive rate (contamination, host artifacts)"
        - "Requires expert interpretation"

      Use_case:
        - "Outbreak investigation"
        - "Unexplained pathology"
        - "Negative screening but clinical suspicion"

      Frequency: "Optional (not routine)"

  #================================
  # STAGE 4: Result Integration
  #================================
  Stage_4_Integration:
    Purpose: "統合判定・定量・品質管理"

    Consolidation:
      Tool: "Custom Python script (pandas + biopython)"

      Input:
        - "Kraken2/Bracken results"
        - "DIAMOND viral results"
        - "MetaPhlAn4 results"
        - "BLAST parasite results"
        - "PERV analysis results"
        - "(Optional) De novo assembly results"

      Process:
        1_Pathogen_list_compilation:
          - "Extract all detected pathogens (>threshold)"
          - "Cross-reference with PMDA 91 list"
          - "Flag PMDA-designated vs. non-designated"

        2_Consensus_calling:
          Logic: "Multi-method agreement"

          High_confidence:
            - "Detected by ≥2 methods"
            - "Read count ≥100"
            - "E-value <1e-10 (BLAST methods)"

          Medium_confidence:
            - "Detected by 1 method"
            - "Read count 50-100"
            - "E-value 1e-5 to 1e-10"

          Low_confidence:
            - "Detected by 1 method"
            - "Read count 10-50"
            - "Requires manual review"

          Discarded:
            - "Read count <10"
            - "Likely contamination/artifact"

        3_Quantification:
          Metric: "Normalized read counts"

          Methods:
            TPM:
              Formula: "(Pathogen_reads / Pathogen_genome_length_kb) / (Total_reads / 1,000,000)"
              Use: "Cross-sample comparison"

            RPM:
              Formula: "(Pathogen_reads / Total_reads) × 1,000,000"
              Use: "Simple abundance"

            Estimated_copies_per_mL:
              Formula: "(Pathogen_reads / Total_reads) × (Input_molecules) / Sample_volume_mL"
              Use: "Absolute quantification (requires spike-in)"

            Spike_in_normalization:
              Spike: "ERCC RNA Spike-In Mix (Ambion)"
              Concentration: "10^4 molecules/μL"
              Volume: "1μL per sample"
              Use: "Accurate absolute quantification"

          Thresholds:
            PMDA_pathogens:
              - Negative: "<10 copies/mL"
              - Positive: "≥10 copies/mL"
              - High_risk: "≥1,000 copies/mL"

            PERV:
              - Acceptable: "<100 copies/mL (DNA level)"
              - Caution: "100-1,000 copies/mL"
              - Alert: ">1,000 copies/mL"

        4_Quality_control:
          Metrics:
            Sequencing_QC:
              - Total_reads: ">20M (PASS/FAIL)"
              - Q30_score: ">85% (PASS/FAIL)"
              - Duplication_rate: "<40% (PASS/WARN/FAIL)"

            Host_depletion_QC:
              - Host_removal: ">95% (PASS/WARN)"
              - Pathogen_enrichment: ">200x (PASS/WARN)"

            Spike_in_QC:
              - ERCC_recovery: "80-120% (PASS)"
              - CV_across_ERCCs: "<20% (PASS)"

            Positive_control:
              - Expected_pathogen_detected: "YES (PASS)"
              - Expected_copy_number: "±20% (PASS)"

            Negative_control:
              - No_pathogen_detected: "YES (PASS)"
              - Contamination_check: "<100 reads total (PASS)"

          Overall_QC_decision:
            - All_PASS: "Release report"
            - Any_FAIL: "Retest / Investigate"
            - WARN_only: "Release with comments"

    Output:
      Files:
        - "pathogen_detection_table.tsv"
        - "quantification_summary.xlsx"
        - "qc_metrics.json"
        - "integrated_report.html"

      Format:
        Table_columns:
          - Pathogen_name
          - PMDA_designated (Y/N)
          - Detection_method(s)
          - Confidence (High/Med/Low)
          - Read_count
          - TPM
          - Estimated_copies_per_mL
          - Taxonomic_lineage
          - GenBank_accession (best match)
          - Alert_level (None/Caution/High)

    Time: "30 min"

  #================================
  # STAGE 5: Reporting & Visualization
  #================================
  Stage_5_Reporting:
    Format: "Multi-format output"

    Report_types:
      1_Executive_summary:
        Format: "PDF (1-2 pages)"
        Language: "Japanese"

        Contents:
          - サンプル情報 (ID, 採取日, etc.)
          - 総合判定: 合格 / 条件付き合格 / 不合格
          - 検出病原体リスト (PMDA指定のみ)
          - リスク評価: 低/中/高
          - 推奨事項
          - 品質管理メトリクス
          - 担当者署名・承認

        Audience: "Clinical decision makers"

      2_Detailed_technical_report:
        Format: "PDF (10-20 pages)"
        Language: "Japanese + English technical terms"

        Contents:
          Section 1: サマリー
          Section 2: 方法論
            - サンプル調製
            - シーケンシング条件
            - 解析パイプライン

          Section 3: 品質管理
            - QCメトリクス (表・グラフ)
            - コントロール結果

          Section 4: 検出結果
            - PMDA 91病原体 スクリーニング結果 (全リスト)
            - 検出病原体詳細 (各病原体ごと)
            - 定量データ (表・グラフ)

          Section 5: PERV特別報告
            - PERV型・サブタイプ
            - ウイルス量
            - 変異解析
            - 系統樹

          Section 6: その他所見
            - 非PMDA病原体検出
            - 未知病原体候補

          Section 7: 考察・推奨

          Appendix:
            - 使用データベースバージョン
            - 参考文献
            - 生データアクセス情報

        Audience: "Scientists, regulators, auditors"

      3_Interactive_HTML_report:
        Tool: "R Markdown + Plotly"

        Features:
          - インタラクティブ図表 (zoom, hover)
          - タブ切り替え (Overview / QC / Pathogens / PERV / etc.)
          - 検索機能
          - データテーブル (sortable, filterable)
          - ヒートマップ (病原体 × サンプル)
          - 系統樹 (interactive phylogenetic tree)

        Audience: "Internal analysis team"

      4_Machine_readable_output:
        Formats:
          - JSON: "Structured data for database"
          - TSV: "Tabular data for Excel"
          - VCF: "Variant calls (PERV)"
          - BIOM: "Microbiome standard format"

        Purpose: "Data archiving, LIMS integration, meta-analysis"

    Visualization:
      Tools:
        - "Python: matplotlib, seaborn, plotly"
        - "R: ggplot2, pheatmap, ggtree"

      Key_figures:
        Fig1_QC_summary:
          - "Read quality distribution (box plot)"
          - "Host depletion efficiency (bar chart)"

        Fig2_Pathogen_overview:
          - "Detected pathogens (bar chart, sorted by abundance)"
          - "PMDA vs non-PMDA (pie chart)"

        Fig3_Quantification:
          - "Pathogen abundance (log scale bar chart)"
          - "Time series (if longitudinal samples)"

        Fig4_PERV_analysis:
          - "PERV type distribution (stacked bar)"
          - "Viral load (scatter plot)"
          - "Phylogenetic tree (ggtree)"
          - "Coverage plot (read depth across genome)"

        Fig5_Taxonomic_composition:
          - "Krona chart (interactive hierarchical)"
          - "Heatmap (pathogen × sample)"

        Fig6_Negative_controls:
          - "Contamination check (bar chart)"

    Automated_delivery:
      Trigger: "Pipeline completion"
      Method: "Email notification + Cloud upload"

      Email:
        To: "Designated recipients"
        Subject: "[PMDA Pathogen Screen] Sample {ID} - {PASS/FAIL}"
        Body:
          - "Executive summary (inline)"
          - "Attachment: PDF report"
          - "Link: Interactive HTML report (cloud)"
          - "Link: Raw data (secure server)"

      Cloud_storage:
        Platform: "AWS S3 / Google Cloud Storage / Azure Blob"
        Retention: "30 years (PMDA requirement)"
        Backup: "3 copies (on-site + 2 cloud regions)"
        Encryption: "AES-256"
        Access_control: "Role-based (RBAC)"

    Time: "15-30 min (automated)"

#================================
# TOTAL PIPELINE TIME
#================================
TOTAL_TIME:
  Compute_time: "8-12 hours"
  Wall_time: "1-2 days (including manual review)"

  Breakdown:
    - QC: "30 min"
    - Host removal: "2 hours"
    - Pathogen detection: "4-6 hours"
    - Integration: "30 min"
    - Reporting: "30 min"
    - Manual review: "2-4 hours"

#================================
# COMPUTATIONAL REQUIREMENTS
#================================
RESOURCES:
  CPU: "64 cores (recommended)"
  RAM: "256GB (recommended, 128GB minimum)"
  Storage:
    - Raw_data: "50-100GB/sample"
    - Intermediate: "100-200GB/sample"
    - Final_output: "10-20GB/sample"
    - Databases: "500GB (total)"

  Network: "High-speed (database downloads)"

  Recommended_setup:
    - "HPC cluster / Cloud (AWS EC2 c6i.16xlarge)"
    - "NAS: 10TB (RAID6)"
    - "SSD: 2TB (temporary files)"

#================================
# COST ESTIMATION (Computational)
#================================
COMPUTATIONAL_COST:
  Phase_1_external:
    - "Included in service fee (no separate charge)"

  Phase_2_shared_HPC:
    - "Cloud compute: 5,000円/sample"
    - "Storage: 500円/sample/month"

  Phase_3_internal:
    - "Electricity: 500円/sample"
    - "Depreciation: 1,000円/sample"
    - "Maintenance: 500円/sample"
    - Total: "2,000円/sample"
```

### 3.2 データベース更新戦略

```yaml
DATABASE_MAINTENANCE:
  Philosophy: "最新の病原体情報を維持しつつ、安定性を確保"

  Update_schedule:
    Kraken2_DB:
      Frequency: "Monthly"
      Source: "NCBI RefSeq (automated rsync)"
      Process:
        - "Download new RefSeq releases"
        - "Filter for bacteria, archaea, viruses, plasmids"
        - "Add PMDA 91 pathogens (if updated)"
        - "Rebuild Kraken2 index"
        - "Validate with test dataset"

      Downtime: "4-6 hours (overnight)"
      Cost: "Computational only"

    RVDB:
      Frequency: "Biannually (Jan, Jul)"
      Source: "https://rvdb.dbi.udel.edu/"
      Process:
        - "Download latest RVDB release"
        - "Merge with custom PERV sequences"
        - "Rebuild DIAMOND index"
        - "Benchmark against previous version"

      Downtime: "2 hours"

    PMDA_Custom_DB:
      Frequency: "Quarterly + ad-hoc"
      Trigger:
        - "New PMDA pathogen designated"
        - "Strain updates (e.g., new influenza)"
        - "PERV variant discovery"

      Process:
        - "Curate sequences from GenBank/DDBJ"
        - "Quality check (completeness, contamination)"
        - "Add to all relevant databases"
        - "Update taxonomy"
        - "Validate detection"

      Downtime: "1 hour"

    Host_genome:
      Frequency: "Annually"
      Trigger:
        - "New pig genome assembly (Ensembl release)"
        - "Modified pig genome updates (CRISPR edits)"

      Process:
        - "Download new assembly"
        - "Incorporate custom modifications (3KO-7TG-59PERV)"
        - "Rebuild BWA/Bowtie2 indices"
        - "Validate host removal efficiency"

      Downtime: "8 hours"

  Version_control:
    System: "Git + Git LFS (Large File Storage)"

    Tracked_files:
      - "Database FASTA files"
      - "Index files (compressed)"
      - "Taxonomy files"
      - "Metadata (versions, sources)"

    Benefits:
      - "Reproducibility (exact database version used)"
      - "Rollback capability (if new DB has issues)"
      - "Audit trail (who updated when)"

    Repository_structure:
      ```
      databases/
      ├── kraken2/
      │   ├── v2023.01/
      │   ├── v2023.02/
      │   └── current -> v2023.02/
      ├── rvdb/
      │   ├── v30.0/
      │   ├── v30.1/
      │   └── current -> v30.1/
      ├── pmda/
      │   ├── v1.0/
      │   ├── v1.1/
      │   └── current -> v1.1/
      └── host/
          ├── sscrofa11.1_base/
          ├── sscrofa11.1_3KO7TG59PERV/
          └── current -> sscrofa11.1_3KO7TG59PERV/
      ```

  Validation_pipeline:
    Purpose: "Ensure new DB doesn't break detection"

    Test_dataset:
      - "Positive controls: 20 samples (known pathogens)"
      - "Negative controls: 5 samples (pathogen-free)"
      - "Spike-in controls: 10 synthetic mixes"

    Metrics:
      - "Sensitivity: >95% (detect known pathogens)"
      - "Specificity: >99% (no false positives)"
      - "Quantification accuracy: ±20%"

    Decision:
      - PASS: "Deploy new DB to production"
      - FAIL: "Investigate issues, rollback if needed"
```

---

## 4. コスト最適化シミュレーション（NGS単独）

### 4.1 Phase別詳細コスト

```yaml
PHASE_1_External_NGS:
  Duration: "Year 0-1"
  Sample_frequency: "6 samples/year"

  Per_sample_cost:
    Sample_prep:
      cfDNA_RNA_extraction: "3,000円 (Zymo kit, internal)"
      Host_depletion_enzymatic: "5,000円 (NEBNext, internal)"
      Host_depletion_hybridization: "8,000円 (QIAseq, internal)"
      Subtotal: "16,000円"

    Library_prep:
      DNA_library: "8,000円 (Illumina DNA Prep, internal)"
      RNA_library: "12,000円 (Illumina RNA Prep, internal)"
      QC_Bioanalyzer: "2,000円"
      Subtotal: "22,000円"

    Sequencing:
      NGS_external: "50,000円 (NovaSeq, 50M reads, outsourced)"
      Data_transfer: "1,000円"
      Subtotal: "51,000円"

    Analysis:
      Bioinformatics: "含まれる (external service)"
      Report_generation: "5,000円 (external)"
      Subtotal: "5,000円"

    QC_Controls:
      Positive_control: "2,000円/run (1 per 6 samples)"
      Negative_control: "500円/run"
      Spike_in: "1,000円"
      Subtotal: "3,500円"

    Personnel:
      Lab_tech: "10,000円/sample (2 hours × 5,000円/hour)"
      Review_time: "5,000円/sample (1 hour × 5,000円/hour)"
      Subtotal: "15,000円"

    Total_per_sample: "112,500円"

  Annual_cost:
    Sample_cost: "112,500円 × 6 = 675,000円"

    Fixed_costs:
      Equipment_depreciation: "500,000円 (extraction equipment)"
      Reagent_stock: "200,000円"
      Training: "300,000円"
      Subtotal: "1,000,000円"

    Total_annual: "1,675,000円"
    Average_per_sample: "279,167円"

  Initial_investment:
    Extraction_equipment: "3,000,000円"
    PCR_equipment: "不要（削除）"
    Analysis_PC: "1,000,000円"
    Training: "500,000円"
    Validation: "3,000,000円"
    Total: "7,500,000円"

PHASE_2_Shared_Sequencing:
  Duration: "Year 1-3"
  Sample_frequency: "12 samples/year"

  Per_sample_cost:
    Sample_prep: "16,000円 (same as Phase 1)"
    Library_prep: "22,000円 (same as Phase 1)"

    Sequencing:
      Shared_facility_fee: "15,000円 (NextSeq, shared use)"
      Reagent_contribution: "6,000円"
      Subtotal: "21,000円"

    Analysis:
      Bioinformatics_internal: "5,000円 (compute + labor)"
      Report_generation: "2,000円 (automated)"
      Subtotal: "7,000円"

    QC_Controls: "3,500円"
    Personnel: "12,000円 (increased efficiency)"

    Total_per_sample: "81,500円"

  Annual_cost:
    Sample_cost: "81,500円 × 12 = 978,000円"

    Fixed_costs:
      Equipment_maintenance: "300,000円"
      Shared_facility_membership: "500,000円/year"
      Bioinformatician_salary: "4,000,000円 (50% allocation)"
      Compute_infrastructure: "1,000,000円 (amortized)"
      Database_maintenance: "200,000円"
      Subtotal: "6,000,000円"

    Total_annual: "6,978,000円"
    Average_per_sample: "581,500円"

  Additional_investment:
    Compute_cluster: "10,000,000円"
    Storage_NAS: "5,000,000円"
    Bioinformatician_hire: "0円 (salary above)"
    Training_advanced: "1,000,000円"
    Total: "16,000,000円"

  ROI:
    Savings_vs_Phase1: "279,167円 - 581,500円 = -302,333円/sample"
    Note: "まだ赤字（投資回収期間）"
    Break_even: "Year 3以降"

PHASE_3_Full_Internal:
  Duration: "Year 3-5"
  Sample_frequency: "24 samples/year"

  Equipment_choice: "MiSeq (cost-optimized)"

  Per_sample_cost:
    Sample_prep: "14,000円 (bulk purchasing discount)"
    Library_prep: "18,000円 (bulk discount)"

    Sequencing:
      MiSeq_reagent: "3,333円 (80,000円/run ÷ 24 samples)"
      Consumables: "1,000円"
      Subtotal: "4,333円"

    Analysis:
      Bioinformatics: "2,000円 (optimized pipeline)"
      Report: "1,000円 (fully automated)"
      Subtotal: "3,000円"

    QC_Controls: "2,000円 (economy of scale)"
    Personnel: "8,000円 (increased efficiency)"

    Total_per_sample: "49,333円"

  Annual_cost:
    Sample_cost: "49,333円 × 24 = 1,184,000円"

    Fixed_costs:
      MiSeq_maintenance: "1,500,000円/year"
      Equipment_depreciation: "1,000,000円/year (MiSeq 15M÷15年)"
      Personnel:
        Lab_tech: "5,000,000円 (full-time)"
        Bioinformatician: "6,000,000円 (75% allocation)"
        Subtotal: "11,000,000円"

      Compute_electricity: "500,000円"
      Database_maintenance: "200,000円"
      Reagent_stock: "500,000円"
      Subtotal: "13,700,000円"

    Total_annual: "14,884,000円"
    Average_per_sample: "620,167円"

  Additional_investment:
    MiSeq: "15,000,000円"
    Automation_equipment: "2,000,000円"
    Total: "17,000,000円"

  Cost_per_sample_marginal: "49,333円"

  ROI:
    Savings_vs_Phase2: "81,500円 - 49,333円 = 32,167円/sample"
    Annual_savings: "32,167円 × 24 = 772,000円"
    Investment_payback: "17,000,000円 ÷ 772,000円 = 22年"

    Note: "検体数24/年では投資回収困難"
    Recommendation: "NextSeq 2000 検討（高throughput）"

PHASE_3_Alternative_NextSeq2000:
  Duration: "Year 3-5"
  Sample_frequency: "24 samples/year → 48 samples/year target"

  Equipment_choice: "NextSeq 2000 (future-proof)"

  Per_sample_cost (48 samples/run):
    Sample_prep: "14,000円"
    Library_prep: "18,000円"

    Sequencing:
      NextSeq2000_reagent: "3,125円 (150,000円/run ÷ 48 samples)"
      Consumables: "1,000円"
      Subtotal: "4,125円"

    Analysis: "3,000円"
    QC_Controls: "1,500円"
    Personnel: "7,000円"

    Total_per_sample: "47,625円"

  Annual_cost (48 samples):
    Sample_cost: "47,625円 × 48 = 2,286,000円"

    Fixed_costs:
      NextSeq2000_maintenance: "2,500,000円/year"
      Equipment_depreciation: "1,667,000円/year (25M÷15年)"
      Personnel: "12,000,000円 (2 staff)"
      Other: "1,500,000円"
      Subtotal: "17,667,000円"

    Total_annual: "19,953,000円"
    Average_per_sample: "415,688円"

  Additional_investment: "25,000,000円"

  ROI:
    Savings_vs_Phase2: "81,500円 - 47,625円 = 33,875円/sample"
    Annual_savings (48 samples): "33,875円 × 48 = 1,626,000円"
    Investment_payback: "25,000,000円 ÷ 1,626,000円 = 15.4年"

    Revised_with_external_revenue:
      External_samples: "24/year (受託サービス)"
      Revenue_per_sample: "150,000円"
      Annual_revenue: "3,600,000円"
      Payback: "25,000,000円 ÷ (1,626,000円 + 3,600,000円) = 4.8年"

    Recommendation: "✓ Phase 4受託サービス前提で投資価値あり"

PHASE_4_Optimized:
  Duration: "Year 5+"
  Sample_frequency: "Internal 36 + External 36 = 72 samples/year"

  Per_sample_cost:
    Sample_prep: "12,000円 (bulk discount)"
    Library_prep: "15,000円 (bulk + automation)"
    Sequencing: "2,083円 (150,000円 ÷ 72 samples)"
    Analysis: "2,000円"
    QC: "1,000円"
    Personnel: "5,000円 (automation)"

    Total_per_sample: "37,083円"

  Annual_cost:
    Sample_cost: "37,083円 × 72 = 2,670,000円"
    Fixed_costs: "18,000,000円"
    Total_annual: "20,670,000円"
    Average_per_sample: "287,083円"

  Revenue_model:
    Internal_samples: "36 × 37,083円 = 1,335,000円 (cost)"
    External_samples: "36 × 150,000円 = 5,400,000円 (revenue)"
    Gross_profit: "5,400,000円 - (36 × 37,083円) - (18,000,000円 × 50%) = 5,400,000円 - 1,335,000円 - 9,000,000円 = -4,935,000円"

    Revised (72 external):
      Revenue: "72 × 150,000円 = 10,800,000円"
      Cost: "72 × 37,083円 + 18,000,000円 = 20,670,000円"
      Profit: "-9,870,000円"

    Note: "外部受託価格を200,000円以上に設定する必要あり"

    Sustainable_model:
      External_price: "200,000円"
      External_samples: "48/year"
      Revenue: "48 × 200,000円 = 9,600,000円"
      Cost: "48 × 37,083円 = 1,780,000円 (marginal)"
      Gross_margin: "9,600,000円 - 1,780,000円 = 7,820,000円"

      Total_revenue: "9,600,000円 + R&D grants (5M) = 14,600,000円"
      Total_cost: "Internal 36 + External 48 = 84 samples"
      Total_cost: "84 × 37,083円 + 18,000,000円 = 21,115,000円"
      Net: "-6,515,000円"

      Conclusion: "Phase 4はR&D投資と位置づけ、外部資金獲得が必須"
```

### 4.2 5年間総コストシミュレーション

```
【シナリオA: 低頻度（年6→12検体）、Phase 2で停止】
┌────────────────────────────────────────┐
│      Year 0   Year 1   Year 2   Year 3   Year 4   Total     │
├────────────────────────────────────────┤
│初期投資  7.5M                                     7.5M      │
│Phase 1   1.7M   1.7M                              3.4M      │
│Phase 2準備      16M                               16M       │
│Phase 2               7.0M    7.0M    7.0M         21M       │
├────────────────────────────────────────┤
│Total     9.2M   17.7M  7.0M    7.0M    7.0M       47.9M     │
│検体数    6      6      12      12      12         48        │
│/検体     1.53M  2.95M  0.58M   0.58M   0.58M      0.998M平均│
└────────────────────────────────────────┘

推奨: 低頻度ならPhase 2維持（外注シーケンス）

【シナリオB: 中頻度（年24検体）、Phase 3 NextSeq移行】
┌────────────────────────────────────────┐
│      Year 0   Year 1   Year 2   Year 3   Year 4   Total     │
├────────────────────────────────────────┤
│初期投資  7.5M                                     7.5M      │
│Phase 1   3.4M                                     3.4M      │
│Phase 2準備 16M                                    16M       │
│Phase 2         7.0M    7.0M                       14M       │
│Phase 3準備              25M                       25M       │
│Phase 3                       20M     20M          40M       │
├────────────────────────────────────────┤
│Total     26.9M 7.0M    32M     20M     20M        105.9M    │
│検体数    24     12      24      24      24         108       │
│/検体     1.12M  0.58M   1.33M   0.83M   0.83M      0.981M平均│
└────────────────────────────────────────┘

推奨: Year 3から大幅コスト削減、Phase 4受託準備

【シナリオC: 高頻度（年36→72検体）、Phase 4展開】
┌────────────────────────────────────────┐
│      Year 0   Year 1   Year 2   Year 3   Year 4   Total     │
├────────────────────────────────────────┤
│初期投資  7.5M                                     7.5M      │
│Phase 1   4.0M                                     4.0M      │
│Phase 2準備 16M                                    16M       │
│Phase 2         7.0M                               7M        │
│Phase 3準備       25M                              25M       │
│Phase 3                20M                         20M       │
│Phase 4準備                   10M                  10M       │
│Phase 4 (internal)                 1.3M   1.3M     2.6M      │
│Phase 4 (external revenue)        -9.6M  -9.6M    -19.2M     │
├────────────────────────────────────────┤
│Total     27.5M 32M     45M     1.7M    -8.3M      98M       │
│Internal検体 36  36     36      36      36         180       │
│/検体     0.76M  0.89M  1.25M   0.05M   -0.23M     0.54M平均 │
└────────────────────────────────────────┘

推奨: 外部受託で収益化、R&D投資として位置づけ
```

### 4.3 Break-Even分析（NGS単独）

```
【内製化判断ポイント】

Phase 1 → Phase 2 移行:
├─ 年間検体数: 6検体以上で有利
├─ 投資額: 16M円
├─ 投資回収: 約2.5年
└─ 推奨: バイオインフォマティクス人材確保できれば即実施

Phase 2 → Phase 3 移行:
├─ 年間検体数: 24検体以上で有利（NextSeq 2000選択時）
├─ 投資額: 25M円
├─ 投資回収:
│  ├─ Internal onlyのみ: 15-20年（非推奨）
│  └─ External受託込み: 5-8年（推奨）
└─ 推奨: 外部受託ビジネスモデル確立後に実施

Phase 3 → Phase 4 移行:
├─ 前提: R&D投資として位置づけ
├─ 必要条件:
│  ├─ 外部受託年間48検体以上確保
│  ├─ 単価200,000円以上設定可能
│  └─ 研究費・補助金獲得（年間5M円以上）
└─ 推奨: 事業計画精査後に慎重判断

【推奨実装パス】
Year 0-1:  Phase 1 (外注NGS)
Year 1-3:  Phase 2 (解析内製化) ← ここまで確実に実施
Year 3-5:  Phase 2継続 or Phase 3検討
           ↓
           受託ビジネス可能性評価
           ↓
          【可能】→ Phase 3 (NextSeq購入)
          【不可】→ Phase 2維持（共同利用継続）
```

---

## 5. 推奨実装戦略（最終版）

### 5.1 推奨アプローチ

```
╔════════════════════════════════════════════════════════════╗
║  NGS Single Platform Strategy (PCR Minimal)                ║
╚════════════════════════════════════════════════════════════╝

【コア戦略】
✓ NGSメタゲノム解析を唯一の検出プラットフォームとする
✓ PCRは原則使用しない（プライマー設計・メンテナンスコスト削減）
✓ 計算科学的バリデーションで精度担保
✓ 段階的内製化によるコスト最適化

【PCR使用が許容される例外ケース】
- 緊急時の迅速スクリーニング（結果まで6時間以内必要）
- NGS装置故障時のバックアップ
- 規制当局が特定病原体でPCR確認を要求する場合のみ

それ以外は全てNGSで完結

【3段階実装】
Phase 1 (Year 0-1): 外注NGS + 内製核酸抽出
├─ 投資: 7.5M円
├─ コスト: 112,500円/検体
├─ 目的: 技術習得、バリデーション
└─ マイルストーン: 解析プログラム取得、10検体経験

Phase 2 (Year 1-3): 解析内製化
├─ 投資: 16M円
├─ コスト: 81,500円 → 27,000円/検体
├─ 目的: データ自主管理、TAT短縮
└─ マイルストーン: バイオインフォマティシャン育成

Phase 2継続 or Phase 3選択 (Year 3-5):
├─ 判断基準:
│  ├─ 年間検体数 ≥24: Phase 3検討
│  ├─ 外部受託見込み: Phase 3推奨
│  └─ 年間検体数 <24: Phase 2維持
│
├─ Phase 3 (装置購入):
│  ├─ 投資: 25M円 (NextSeq 2000)
│  ├─ コスト: 47,625円/検体（48検体時）
│  └─ ROI: 外部受託込みで5年回収
│
└─ Phase 2継続:
   └─ 共同利用施設を継続使用（柔軟性高い）
```

### 5.2 即座のアクション（優先順位付き）

```
【Week 1: プロジェクト始動】
Priority 1 (Critical):
□ プロジェクトチーム編成
  └─ PM, ラボ責任者, バイオインフォ担当（外部コンサル可）
□ 予算承認プロセス開始
  └─ Phase 1: 7.5M円（初期投資）
  └─ Phase 1: 1.7M円（年間運営、6検体）
□ 外注先NGSサービス 3社に見積依頼
  候補:
  ├─ タカラバイオ（NGS受託No.1、国内、信頼性高）
  ├─ Macrogen Japan（価格競争力、韓国本社、実績豊富）
  └─ Azenta Life Sciences（旧GENEWIZ、米系、先進的）

Priority 2 (High):
□ 装置選定開始（核酸抽出、QC）
  ├─ 自動核酸抽出装置: Zymo QuickDNA/RNA 96 Kit対応
  └─ QC: Agilent TapeStation 4150 or 4200
□ PMDA事前相談予約
  └─ バリデーション計画の事前確認

【Month 1: セットアップ】
Priority 1:
□ 外注先選定・契約締結
  └─ 解析プログラムコード取得条項を明記
□ 装置発注
  ├─ 核酸抽出装置: 3M円
  ├─ TapeStation: 2M円
  └─ 分析PC (ハイスペック): 1M円
□ ラボスペース準備
  └─ -80°C冷凍庫、クリーンベンチ確認

Priority 2:
□ SOP作成開始
  ├─ cfDNA/RNA抽出プロトコル
  ├─ サンプル管理
  └─ データ記録

【Month 2-3: トレーニング & パイロット】
□ ラボスタッフトレーニング
  ├─ cfDNA/RNA抽出技術
  ├─ ライブラリQC
  └─ データ解釈基礎
□ パイロット実験 (n=3サンプル)
  ├─ 陽性コントロール（既知病原体スパイク）
  ├─ 陰性コントロール（SPFブタ血漿）
  └─ 実サンプル（候補ドナーブタ）
□ 外注先から解析コード受領
  └─ ローカル環境で再現テスト

【Month 4-6: バリデーション準備】
□ 標準サンプル調製（LOD決定用）
  └─ 主要10病原体 × 6濃度 = 60サンプル
□ バリデーション実験計画書作成
  └─ PMDAガイドライン準拠
□ 倫理審査・動物実験承認
  └─ サンプル採取プロトコル

【Month 7-12: 分析的バリデーション】
□ LOD決定実験
□ 再現性検証（operator間、日間、装置間）
□ 特異性試験（交差反応性）
□ データ解析SOPファイナライズ
□ バリデーションレポート作成
  └─ PMDA提出用

【Month 13-: Phase 2準備 & 移行判断】
□ Phase 1レビュー（Go/No-Go判断）
  ├─ 技術習熟度評価
  ├─ コスト実績確認
  └─ Phase 2投資対効果再計算
□ バイオインフォマティシャン採用開始
  ├─ スキル要件: NGS解析、Python/R、Linux
  └─ 給与レンジ: 5-7M円/年
□ 計算インフラ設計
  └─ クラウド vs オンプレミス検討
```

### 5.3 成功のための重要ポイント（NGS特化版）

```
【Critical Success Factors - NGS Single Platform】

1. データベース品質の維持
   ├─ PMDA 91病原体の網羅的カバー
   ├─ 定期更新（月次〜四半期）
   ├─ バージョン管理の徹底
   └─ バリデーション再実施（DB更新時）

   Risk: データベース古い → 新規株・変異株の見逃し
   Mitigation: 自動更新スクリプト + 検証パイプライン

2. シーケンス深度の確保
   ├─ 最低50M reads/sample (Phase 1-2)
   ├─ ホスト除去効率 >95%
   ├─ 病原体エンリッチメント >200x
   └─ LOD達成: 10 copies/mL

   Risk: リード数不足 → 偽陰性
   Mitigation: QC gate設定、不合格サンプルの再シーケンス

3. バイオインフォマティクス人材
   ├─ 専任1名以上（Phase 2-）
   ├─ スキル: NGS, Linux, Python/R, 統計学
   ├─ 継続的教育・学会参加
   └─ 外部専門家との連携

   Risk: 人材流出、技術陳腐化
   Mitigation: 複数スタッフのクロストレーニング、SOPの充実化

4. 計算インフラの信頼性
   ├─ 24/7稼働可能なサーバー
   ├─ 自動バックアップ（3拠点）
   ├─ 災害復旧計画（DR）
   └─ セキュリティ（データ暗号化、アクセス制御）

   Risk: システム障害 → 解析遅延
   Mitigation: 冗長構成、クラウドバックアップ

5. PCRフリー戦略の堅持
   ├─ NGS結果の信頼性を計算科学的に担保
   ├─ 多重検証（Kraken2 + BLAST + MetaPhlAn）
   ├─ 偽陽性フィルタリング（リード数、E-value）
   └─ 陽性結果の手動レビュー

   Risk: PCR確認なしへの不安 → 採用躊躇
   Mitigation: 徹底したバリデーション、論文発表で信頼性証明

6. 規制対応の継続的モニタリング
   ├─ PMDAガイドライン変更を追跡
   ├─ FDA/EMA等の国際動向も参照
   ├─ 年1回のレビューミーティング
   └─ 必要に応じてバリデーション追加

   Risk: 規制変更 → システム不適合
   Mitigation: 柔軟なシステム設計、余裕のあるバリデーション

7. 段階的投資の厳格な判断
   ├─ 各フェーズ終了時にGo/No-Go判断
   ├─ KPI達成度評価（技術、コスト、品質）
   ├─ 次フェーズのROI再計算
   └─ 外部環境変化の考慮（装置価格、受託相場等）

   Risk: 過剰投資、投資タイミング誤り
   Mitigation: マイルストーンゲート、外部専門家レビュー
```

---

## 6. 結論

### 6.1 最終推奨事項

```
╔══════════════════════════════════════════════════════════════╗
║  PCR Minimal, NGS-Centric Strategy for PMDA 91 Pathogens    ║
╚══════════════════════════════════════════════════════════════╝

【戦略サマリー】
- NGSメタゲノム解析を唯一の主軸プラットフォームとする
- PCRは使用しない（例外: 緊急時、装置故障時のみ）
- 段階的内製化により長期的コスト削減

【コスト予測】
Phase 1 (Year 0-1):   112,500円/検体
Phase 2 (Year 1-3):    27,000円/検体  (-76% vs Phase 1)
Phase 3 (Year 3-5):    47,625円/検体  (48検体/年時)
Long-term:             37,000円/検体  (最適化後)

vs. PCR併用戦略 v1.0:  180,000円/検体 → 50,000円/検体
差額: NGS単独で30-40%さらに削減

【技術的優位性】
✓ 91種類すべて + 未知病原体を同時検出
✓ 変異株・新規株に自動対応（プライマー再設計不要）
✓ 系統解析・定量を統合的に実施
✓ データの二次利用価値（研究、AI、疫学）

【リスク】
△ 初期投資大（7.5M → 16M → 25M円）
△ バイオインフォマティクス専門人材必須
△ NGS装置故障時のバックアップ必要

【推奨実装パス】
┌─────────────────────────────────┐
│ Year 0-1:  Phase 1 (外注NGS)               │
│            ├─ 技術習得、バリデーション      │
│            └─ 解析プログラム取得            │
│                                             │
│ Year 1-3:  Phase 2 (解析内製化)            │
│            ├─ バイオインフォ人材確保       │
│            ├─ 計算インフラ構築             │
│            └─ コスト76%削減達成            │
│                                             │
│ Year 3-5:  Phase 2継続 or Phase 3移行      │
│            Decision point:                  │
│            ├─ 年間24検体以上: Phase 3      │
│            ├─ 外部受託可能: Phase 3        │
│            └─ それ以外: Phase 2継続        │
│                                             │
│ Year 5+:   Phase 4 (Optional)              │
│            ├─ R&D投資として位置づけ        │
│            ├─ 外部受託ビジネス化           │
│            └─ AI/機械学習統合              │
└─────────────────────────────────┘

【投資意思決定】
Approve Phase 1: ✓ 推奨（必須投資）
  └─ 7.5M円、確実なROI

Approve Phase 2: ✓ 強く推奨
  └─ 16M円、2.5年で回収、データ主権確保

Approve Phase 3: △ 条件付き推奨
  └─ 25M円、検体数・受託ビジネス次第
  └─ Phase 2終了時に再評価

Approve Phase 4: ○ オプション
  └─ 20M円、戦略的R&D投資、外部資金前提
```

### 6.2 Next Steps (今週実施)

```
【即座のアクション】
□ 本提案書を経営層に提出（意思決定依頼）
□ Phase 1予算承認プロセス開始（7.5M + 1.7M円）
□ プロジェクトマネージャー任命
□ 外注NGS企業3社に見積依頼（RFP送付）
□ PMDAとの事前相談予約

【承認後1週間以内】
□ キックオフミーティング開催
□ 装置選定・発注
□ ラボスペース確保
□ SOPドラフト作成開始

【承認後1ヶ月以内】
□ 外注先選定・契約締結
□ トレーニングプログラム開始
□ パイロット実験計画finalize
```

---

**Document Version**: 2.0 (PCR Minimal Strategy)
**Created**: 2025-10-08
**Author**: Pathogen Detection System Design Team
**Next Review**: Project Month 3, 6, 12
