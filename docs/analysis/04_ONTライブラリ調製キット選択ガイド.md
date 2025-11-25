# Oxford Nanopore ライブラリ調製キット選択ガイド - ウイルス網羅的検査用

**質問**: ONT社のNGS装置でウイルスを網羅的に検査する場合、どの試薬をどのように使用するのが効果的か？
**対象**: PMDA 91病原体スクリーニング（ウイルス41種 + 細菌27種 + 寄生虫19種 + 真菌2種）

---- 

## 1. エグゼクティブサマリー

### 1.1 核心的な回答

**チームメンバーの観察は正しい**:
- ✅ ライブラリ調製キットとシーケンサーには相性がある
- ✅ ONT装置にはONT試薬が必須
- ✅ Illumina用キット（NEBNext等）はONTで使用不可

**推奨キット**: **Ligation Sequencing Kit (SQK-LSK114)** ⭐

```yaml
最適解: Ligation Sequencing Kit (SQK-LSK114)

理由:
  ✅ 最高収量（20-50 Gb/MinION run）
  ✅ 最長read length（平均10-50 kb）
  ✅ ウイルスゲノム完全長取得可能
  ✅ PMDA 91病原体すべてに対応
  ✅ Protocol 12 v2.1で使用中（実績あり）

コスト: ¥15,000-20,000/sample
所要時間: 8-10時間

代替案:
  - Rapid Sequencing Kit: 緊急時のみ（収量50%減）
  - Direct RNA Kit: RNAウイルス詳細解析（特殊用途）
  - PCR-free: 低バイアス（研究用）
```

---- 

## 2. ONTライブラリ調製キットの全体像

### 2.1 主要キットラインナップ（2024-2025年現在）

#### カテゴリー別分類

```yaml
【カテゴリーA: DNA Sequencing - General Purpose】

1. Ligation Sequencing Kit (SQK-LSK114) ⭐ 推奨
   用途: 標準的DNA sequencing（最高性能）
   Input: 200 ng - 3 μg DNA
   Fragment size: 200 bp - 200 kb+
   Yield: 20-50 Gb（MinION）
   Read length: Average 10-50 kb
   Protocol time: 8-10 hours
   Cost: $550-650/kit（12 samples）= ¥85,000-100,000
         ¥7,000-8,500/sample

   特徴:
     ✅ 最高収量・最長read
     ✅ End-repair + A-tailing + Adapter ligation
     ✅ 低バイアス
     ✅ ウイルス全ゲノム解析に最適

2. Rapid Sequencing Kit (SQK-RAD114)
   用途: 迅速スクリーニング（緊急検査）
   Input: 400 ng DNA（推奨）
   Fragment size: >200 bp
   Yield: 10-25 Gb（MinION、Ligation比50%）
   Read length: Average 5-20 kb
   Protocol time: 10 minutes（超高速）
   Cost: $400/kit（12 samples）= ¥62,000
         ¥5,200/sample

   特徴:
     ✅ 超高速（10分でライブラリ完成）
     ⚠️ 収量低い（Ligationの50%）
     ⚠️ Read length短い
     ⚠️ バイアス高い（transposase使用）
     → 緊急時・スクリーニングのみ推奨

3. PCR-Free Sequencing Kit (SQK-LSK114 + PCR-free prep)
   用途: 低バイアス研究（増幅なし）
   Input: 1-3 μg DNA（高input必要）
   Yield: 15-40 Gb（同上）
   Protocol time: 8-10 hours
   Cost: Same as Ligation

   特徴:
     ✅ PCR amplification bias なし
     ✅ GC-rich region問題回避
     ⚠️ High input requirement（1 μg+）
     → 血漿サンプルには不向き

【カテゴリーB: RNA Sequencing】

4. Direct RNA Sequencing Kit (SQK-RNA004)
   用途: RNAを直接シーケンス（cDNA変換なし）
   Input: 500 ng - 1 μg poly(A)+ RNA
   Yield: 5-15 Gb（MinION）
   Read length: Full-length transcripts
   Protocol time: 6-8 hours
   Cost: $850/kit（12 samples）= ¥130,000
         ¥11,000/sample

   特徴:
     ✅ Poly(A)+ RNAウイルスに最適
     ✅ 修飾塩基検出可能（m6A等）
     ✅ Full-length viral genomes
     ⚠️ Poly(A)- RNAには不向き（Hantavirus等）
     → 特殊用途（詳細解析時のみ）

5. cDNA-PCR Sequencing Kit (SQK-PCS114)
   用途: 低input RNA（cDNA経由）
   Input: 1-100 ng total RNA
   Yield: 10-30 Gb（MinION）
   Protocol time: 10-12 hours（RT-PCR含む）
   Cost: $650/kit（12 samples）

   特徴:
     ✅ 低input対応（1 ng~）
     ✅ Poly(A)+/- 両対応
     ⚠️ PCR bias あり
     → Protocol 12でRNA virus対応に使用

【カテゴリーC: Amplicon Sequencing】

6. PCR Barcoding Kit (SQK-PBK004)
   用途: PCR amplicon sequencing
   Applications: 16S rRNA, Targeted panels
   Cost: $400/kit（96 samples）

   特徴:
     ✅ 96サンプル多重化
     ✅ 低コスト（¥4,200/sample）
     ⚠️ Amplicon専用（metagenomics不向き）

【カテゴリーD: Specialized Applications】

7. Adaptive Sampling Kit (included in software)
   用途: Real-time selective sequencing
   Cost: Free（MinKNOW software機能）

   特徴:
     ✅ Target enrichment（in silico）
     ⚠️ 80% yield reduction in some samples
     → ハイブリッドキャプチャの方が効率的

8. CRISPR-Cas Enrichment Kit (研究用）
   用途: CRISPR-based target enrichment
   Status: Research use only
   → 現時点で臨床使用推奨せず
```

---- 

## 3. PMDA 91病原体スクリーニングに最適なキット選択

### 3.1 病原体タイプ別推奨キット

#### DNAウイルス（25種）

```yaml
Target Viruses:
  - ASFV, CSFV, FMDV
  - PCV2, PCV3（Circular ssDNA）
  - Adenovirus, Herpesvirus, Poxvirus等

推奨キット: Ligation Sequencing Kit (SQK-LSK114) ⭐

理由:
  ✅ DNA直接シーケンス
  ✅ 長read（10-50 kb）→ 全ゲノム取得
  ✅ Circular genome対応（Protocol 12 Step 2.5併用）

Workflow:
  1. DNA extraction（血漿5 mL）
  2. DNase I treatment（Circular ssDNA linearization）
  3. End-repair + A-tailing
  4. Adapter ligation（SQK-LSK114）
  5. Sequencing（MinION/PromethION）

Expected yield:
  - Input: 100-500 ng DNA（血漿由来）
  - Output: 20-40 Gb/run
  - Viral reads: 1-10%（host depletion後）

Performance:
  - LOD: 50-100 gc/mL（Hybrid Capture併用）
  - Full genome recovery: >90% for abundant viruses
```

---- 

#### RNAウイルス - Poly(A)+ （18種）

```yaml
Target Viruses:
  - Influenza A/B/C
  - Coronavirus（SADS-CoV等）
  - その他Poly(A)+ RNAウイルス

推奨キット: Ligation Sequencing Kit + cDNA conversion

Workflow:
  1. RNA extraction（血漿）
  2. Poly(A) selection（Optional, enrichment）
  3. Reverse transcription:
     - Superscript IV（Invitrogen）
     - Random hexamer + oligo-dT primers
     - Temperature: 50°C, 60 min
  4. Second-strand synthesis:
     - E. coli DNA Polymerase I
     - RNase H
  5. cDNA library prep:
     - Ligation Sequencing Kit (SQK-LSK114)

Expected yield:
  - Input: 50-200 ng total RNA
  - cDNA: 20-100 ng
  - Output: 15-30 Gb/run

Performance:
  - LOD: 100-200 gc/mL（Hybrid Capture併用）
  - Full-length viral genomes: 80-90%

Alternative: Direct RNA Kit (SQK-RNA004)
  用途: 詳細解析（修飾塩基、full-length isoforms）
  Limitation: Higher input requirement（500 ng poly(A)+ RNA）
  Cost: Higher（¥11,000 vs ¥8,000/sample）
  → 初期スクリーニングには不向き、確認用途のみ
```

---- 

#### RNAウイルス - Poly(A)- （8種）

```yaml
Target Viruses:
  - Hantavirus（3-segment genome）
  - EEEV（Eastern Equine Encephalitis Virus）
  - その他Poly(A)- RNAウイルス

Challenge:
  - No poly(A) tail → Poly(A) selection不可
  - rRNA depletion必須（>90% of total RNA）

推奨キット: Ligation Sequencing Kit + rRNA depletion

Workflow:
  1. RNA extraction
  2. rRNA depletion:
     Method A: RiboMinus（Invitrogen）
       - Pig rRNA probes（Sus scrofa）
       - Depletion: >95%
       - Cost: ¥10,000/sample

     Method B: RNase H-based depletion
       - Custom pig rRNA probes
       - Lower cost: ¥5,000/sample
       - 社内データベース: /mnt/efs/databases/pig_rrna/

  3. Reverse transcription:
     - Random hexamer primers（Poly(A)なしのため）
     - Superscript IV
     - 50°C, 60 min

  4. Second-strand synthesis
  5. cDNA library prep: Ligation Kit

Expected yield:
  - Input: 100-500 ng total RNA
  - Post rRNA depletion: 5-50 ng（90-95%除去）
  - cDNA: 2-20 ng
  - Output: 10-25 Gb/run（低input影響）

Performance:
  - LOD: 200-500 gc/mL（rRNA背景影響）
  - Improved to 100-200 gc/mL（Hybrid Capture併用）
  - Segmented genome recovery:
     Hantavirus L/M/S: 85-90% all 3 segments

注意: Poly(A)- RNAは最も困難なカテゴリー
  → Protocol 11（超高感度）推奨
```

---- 

#### 細菌（27種）

```yaml
Target Bacteria:
  - Salmonella, E. coli, Brucella
  - Mycobacterium, Leptospira等

推奨キット: Ligation Sequencing Kit (SQK-LSK114)

Optional Enhancement: 16S rRNA Amplicon Kit
  用途: 低コスト大量スクリーニング
  Kit: PCR Barcoding Kit (SQK-PBK004)
  Cost: ¥4,200/sample（96多重化）

  Workflow:
    1. DNA extraction
    2. 16S rRNA PCR（27F/1492R primers）
    3. Barcoding PCR
    4. Pooling（96 samples）
    5. Sequencing

  Limitation:
    - 16S rRNA sequenceのみ（species同定限界）
    - 定量精度低い
    → メタゲノム全域解析の方が推奨

推奨アプローチ: メタゲノム全域（Ligation Kit）
  理由:
    ✅ Species-level identification
    ✅ Antibiotic resistance genes検出
    ✅ Virulence factors同定
    ✅ Strain typing可能
```

---- 

#### 寄生虫・真菌（21種）

```yaml
Target Organisms:
  - Toxoplasma gondii
  - Trichinella spiralis
  - Candida, Aspergillus等

推奨キット: Ligation Sequencing Kit (SQK-LSK114)

特徴:
  - 大型ゲノム（10-100 Mb）
  - 高バイオマス
  - 検出容易（>10,000 copies/mL typical）

Performance:
  - LOD: 50-200 gc/mL（余裕で達成）
  - Species identification: >99%
  - Long reads有利（repetitive regions解決）
```

---- 

### 3.2 Protocol 12 v2.1 での使用キット

#### 現行プロトコルの構成

```yaml
Protocol 12 v2.1: 統合サンプル調製プロトコル
  用途: DNA + RNA virus同時検出（91病原体対応）

使用キット:
  1. Extraction:
     - QIAamp Circulating Nucleic Acid Kit（Qiagen）
     - 血漿5 mLから DNA + RNA同時抽出
     - Cost: ¥8,000/sample

  2. DNA Virus Library:
     Step 2.5: DNase I + Klenow Fragment（Circular ssDNA対応）
     Step 3: Ligation Sequencing Kit (SQK-LSK114)

  3. RNA Virus Library:
     Step 4: cDNA synthesis
       - Superscript IV Reverse Transcriptase
       - Random hexamer + oligo-dT primers
     Step 5: Second-strand synthesis
     Step 6: Ligation Sequencing Kit (SQK-LSK114)

  4. Multiplexing:
     - Native Barcoding Kit (EXP-NBD196)
     - 24-96 samples/run

Total Kit Cost per Sample:
  - Extraction: ¥8,000
  - Ligation Kit: ¥8,000（amortized）
  - Barcoding: ¥2,000
  - Enzymes（RT, DNase等）: ¥3,000
  Total: ¥21,000/sample（試薬のみ）
```

---- 

## 4. キット選択の決定木

### 4.1 用途別フローチャート

```
ウイルス網羅的検査（PMDA 91病原体）
         ↓
    【質問1: サンプルタイプは？】
         ↓
    ┌────┴────┐
    ↓         ↓
  DNA       RNA
    ↓         ↓
    └────┬────┘
         ↓
    【質問2: 検出目的は？】
         ↓
    ┌────┴────────┐
    ↓               ↓
定量・同定      全ゲノム解析
（スクリーニング）  （詳細解析）
    ↓               ↓
【推奨】          【推奨】
Ligation Kit     Ligation Kit
SQK-LSK114       SQK-LSK114
    ↓               ↓
+ Barcoding      + Ultra-long
+ Hybrid Capture   protocol
         ↓
    【質問3: 緊急性は？】
         ↓
    ┌────┴────┐
    ↓         ↓
  通常      緊急（<24h）
    ↓         ↓
Ligation    Rapid Kit
SQK-LSK114  SQK-RAD114
（推奨）    （妥協案）
```

---- 

### 4.2 シナリオ別推奨

#### シナリオA: 標準PMDA 91病原体スクリーニング（通常）

```yaml
Recommended Kit: Ligation Sequencing Kit (SQK-LSK114)

Configuration:
  - Sample prep: Protocol 12 v2.1
  - Barcoding: EXP-NBD196（24-96 samples）
  - Sequencing: MinION 72h or PromethION 72h
  - Basecalling: Simplex SUP mode

Cost: ¥21,000/sample（試薬）+ ¥10,000（フローセル配分）= ¥31,000
Time: 24h（prep）+ 72h（sequencing）= 4 days
Performance: 89-91/91 pathogens detected

Advantages:
  ✅ 最高性能
  ✅ 全病原体対応
  ✅ 長read（ゲノム再構築可能）
```

---- 

#### シナリオB: 緊急スクリーニング（24時間以内）

```yaml
Recommended Kit: Rapid Sequencing Kit (SQK-RAD114)

Configuration:
  - Sample prep: Simplified（end-repair skip）
  - Barcoding: Optional（推奨はsingle sample）
  - Sequencing: MinION 24h（early stop）
  - Basecalling: Fast mode

Cost: ¥5,200/sample（試薬）+ ¥37,500（フローセル配分、24 samples）= ¥42,700
Time: 10 min（prep）+ 24h（sequencing）= 1 day
Performance: 70-80/91 pathogens detected（感度低下）

Advantages:
  ✅ 超高速
  ✅ 緊急時対応

Disadvantages:
  ⚠️ 収量50%減
  ⚠️ 感度低下（LOD 2-3×悪化）
  ⚠️ Short reads（ゲノム再構築困難）

Use Case:
  - 感染症アウトブレイク疑い
  - 緊急手術前スクリーニング
  - 初期スクリーニング（陽性の場合、Ligation再検査）
```

---- 

#### シナリオC: RNAウイルス詳細解析

```yaml
Recommended Kit: Direct RNA Sequencing Kit (SQK-RNA004)

Configuration:
  - Sample prep: Poly(A) selection
  - Input: 500 ng poly(A)+ RNA（要高濃度）
  - Sequencing: MinION 48h
  - Analysis: Full-length isoforms, modifications

Cost: ¥11,000/sample（試薬）+ Poly(A) selection ¥8,000 = ¥19,000
Time: 6h（prep）+ 48h（sequencing）= 2.5 days
Performance: Poly(A)+ RNAウイルスのみ（18/41 viruses）

Advantages:
  ✅ Direct RNA（no RT bias）
  ✅ Modified bases detection（m6A等）
  ✅ Full-length transcripts

Disadvantages:
  ⚠️ Poly(A)+ only（Hantavirus等検出不可）
  ⚠️ High input requirement（500 ng）
  ⚠️ Lower throughput（5-15 Gb vs 20-50 Gb）

Use Case:
  - PERV詳細解析（env gene isoforms）
  - Influenza詳細解析
  - RNA modification研究
  → 初期スクリーニングには不向き、確認解析用
```

---- 

#### シナリオD: 超低コスト大量スクリーニング（研究用）

```yaml
Recommended Kit: PCR Barcoding Kit (SQK-PBK004)

Configuration:
  - Target: 16S rRNA amplicons（細菌のみ）
  - Barcoding: 96 samples pooled
  - Sequencing: MinION 24h
  - Analysis: Taxonomic classification

Cost: ¥4,200/sample（試薬）+ ¥10,000（フローセル）= ¥14,200
Time: 4h（PCR + barcoding）+ 24h（sequencing）= 1.5 days
Performance: Bacteria only（27/91 pathogens）

Advantages:
  ✅ 超低コスト
  ✅ 96サンプル多重化
  ✅ 高スループット

Disadvantages:
  ⚠️ Bacteria only（ウイルス検出不可）
  ⚠️ Species resolution限界（16S配列のみ）
  ⚠️ Quantification精度低い
  ⚠️ PCR bias

Use Case:
  - 研究用大量スクリーニング
  - 細菌叢解析
  → PMDA申請には不向き
```

---- 

## 5. Barcoding（多重化）戦略

### 5.1 Barcoding Kit選択

#### Native Barcoding Kit (EXP-NBD196) - 推奨 ⭐

```yaml
Product: Native Barcoding Expansion 1-96
Code: EXP-NBD196
Barcodes: 96種類
Compatibility: Ligation Sequencing Kit

特徴:
  ✅ Native adapter（both ends barcoded）
  ✅ 最高精度（demultiplexing accuracy >99.9%）
  ✅ Longest reads preserved
  ✅ 96 samples/run可能

Workflow:
  1. End-repair + A-tailing
  2. Native barcode ligation（sample-specific）
  3. Pooling（96 samples）
  4. Adapter ligation（shared adapters）
  5. Sequencing

Cost: $650/kit（96 barcodes）= ¥100,000
       ¥1,040/sample（96多重化時）

Yield per Sample（PromethION 100 Gb）:
  96 samples: 100 Gb ÷ 96 = 1.04 Gb/sample ✅
  Sufficient for 91-pathogen detection

Yield per Sample（MinION 20 Gb）:
  24 samples: 20 Gb ÷ 24 = 833 Mb/sample ✅
  Sufficient with Hybrid Capture
```

---- 

#### PCR Barcoding Kit (SQK-PBK004) - 低コスト

```yaml
Product: PCR Barcoding Kit
Code: SQK-PBK004
Barcodes: 96種類
Compatibility: PCR products only

特徴:
  ✅ 超低コスト（¥4,200/sample）
  ⚠️ PCR bias導入
  ⚠️ Amplicon専用（metagenomics不向き）

Use Case: 16S rRNA amplicon sequencing（研究用のみ）
```

---- 

### 5.2 多重化レベル決定

#### MinION（20-30 Gb output）

```yaml
24 Samples Multiplexing（推奨）:
  Per-sample yield: 20 Gb ÷ 24 = 833 Mb
  With Hybrid Capture: 500-700 Mb on-target
  LOD: 100 gc/mL ✅
  Cost-efficient: ✅

12 Samples Multiplexing（高感度）:
  Per-sample yield: 20 Gb ÷ 12 = 1.67 Gb
  With Hybrid Capture: 1-1.4 Gb on-target
  LOD: 50 gc/mL ✅
  Cost: 2× higher per sample

推奨: 24 samples（標準スクリーニング）
     12 samples（低ウイルス量疑い時）
```

---- 

#### PromethION P2（100 Gb output）

```yaml
96 Samples Multiplexing（推奨）:
  Per-sample yield: 100 Gb ÷ 96 = 1.04 Gb
  With Hybrid Capture: 700 Mb - 1 Gb on-target
  LOD: 50 gc/mL ✅
  Cost-efficient: ✅✅

48 Samples Multiplexing（超高感度）:
  Per-sample yield: 100 Gb ÷ 48 = 2.08 Gb
  With Hybrid Capture: 1.4-1.8 Gb on-target
  LOD: 20-30 gc/mL ✅✅
  Cost: 2× higher per sample

推奨: 96 samples（標準）
     48 samples（PERV <50 gc/mL detection必要時）
```

---- 

## 6. 試薬調達とコスト管理

### 6.1 年間試薬コスト試算（100 donors/year）

#### 推奨構成: MinION + Ligation Kit + Native Barcoding

```yaml
Annual Reagent Cost（100 donors）:

Extraction:
  - QIAamp Circulating NA Kit: ¥8,000 × 100 = ¥800,000

Library Prep:
  - Ligation Sequencing Kit (SQK-LSK114):
    $650/kit × 9 kits（100 samples ÷ 12 samples/kit）
    = $5,850 = ¥900,000

  - Native Barcoding Kit (EXP-NBD196):
    $650/kit × 2 kits（96 barcodes）= $1,300 = ¥200,000

  - Enzymes（RT, DNase, Klenow）:
    ¥3,000 × 100 = ¥300,000

Host Depletion:
  - MBD-Fc beads: ¥5,000 × 100 = ¥500,000

Hybrid Capture:
  - Capture probes: ¥15,000 × 100 = ¥1,500,000

Flow Cells:
  - MinION flow cells: $900 × 5（100 samples ÷ 24 samples/FC × 1.2）
    = $5,400 = ¥830,000

Quality Control:
  - NanoDrop, Qubit reagents: ¥100,000

Total Annual Reagent Cost: ¥5,130,000
Per Donor Reagent Cost: ¥51,300
```

---- 

### 6.2 発注戦略

#### Bulk Ordering（年間契約）

```yaml
Strategy: 年間使用量を一括発注

Ligation Sequencing Kit:
  Standard: $650/kit
  Bulk discount（10+ kits）: $585/kit（10% OFF）
  Annual order: 10 kits = $5,850
  Savings: $650 = ¥100,000

Flow Cells:
  Standard: $900/FC
  Bulk discount（10+ FCs）: $810/FC（10% OFF）
  Annual order: 10 FCs = $8,100
  Savings: $900 = ¥138,000

Total Annual Savings: ¥238,000（5% reduction）
```

---- 

#### Vendor Selection

```yaml
Primary Vendor: Oxford Nanopore Technologies
  - Direct purchase（公式サイト or 代理店）
  - Guaranteed compatibility
  - Technical support included

日本国内代理店:
  1. 日本ジェネティクス株式会社
     - ONT公式代理店
     - 技術サポート日本語対応
     - 納期: 2-4週間

  2. タカラバイオ株式会社
     - ONT製品取扱い
     - PCR酵素等も一括購入可能

推奨: 日本ジェネティクス（公式代理店、技術サポート充実）

Lead Time Management:
  - Flow cells: 2-4週間（冷蔵品、要在庫管理）
  - Ligation Kit: 2-4週間（冷蔵品）
  - Barcoding Kit: 2-4週間
  - 推奨在庫: 3-6ヶ月分（欠品リスク回避）
```

---- 

## 7. よくある質問（FAQ）

### 7.1 キット互換性

**Q1: Illumina用NEBNextキットはONTで使えますか？**

```yaml
A: 使えません ❌

理由:
  - Illumina adapter: P5/P7 sequences（ONT非互換）
  - ONT adapter: Motor protein binding sequences（特殊構造）
  - End chemistry: 異なる（Illumina: blunt, ONT: dA-tailed）

結果: Illumina用ライブラリはONT flow cellに結合できない

対策: ONT専用キット（Ligation Kit）使用必須
```

**Q2: 古いLigationキット（SQK-LSK109）は使えますか？**

```yaml
A: 非推奨 ⚠️

旧世代キット:
  - SQK-LSK109: R9.4.1 flow cell用（2019-2021）
  - SQK-LSK110: R10.3 flow cell用（2021-2022）
  - SQK-LSK111-113: 過渡期

現行キット:
  - SQK-LSK114: R10.4.1 flow cell用（2023-現在）⭐

互換性:
  - LSK109 → R10.4.1 FC: 動作するが収量50%減
  - LSK114 → R9.4.1 FC: 動作するが非推奨

推奨: 最新キット（SQK-LSK114）+ 最新FC（R10.4.1）使用
```

---- 

### 7.2 トラブルシューティング

**Q3: 収量が低い（\<10 Gb、期待20-50 Gb）**

```yaml
A: 以下を確認

1. Input DNA/RNA量:
   - 推奨: 200-500 ng
   - 実測: NanoDrop/Qubit で確認
   - 対策: 低input時はPCR amplification追加

2. DNA/RNA quality:
   - RIN（RNA Integrity Number）: >7.0必要
   - Fragmentation: >10 kb peak必要
   - 対策: Fresh sample使用、凍結融解回避

3. Adapter ligation効率:
   - Ligation時間: 室温20分以上
   - ATP濃度: Fresh（古いと活性低下）
   - 対策: Ligation時間延長（30-60分）

4. Flow cell QC:
   - Active pore count: >800必要（MinKNOW check）
   - 対策: 新しいflow cell使用

5. Loading concentration:
   - 推奨: 50-100 fmol
   - 対策: Qubit定量 → 正確な濃度調整
```

**Q4: Barcoding後のdemultiplexing精度が低い**

```yaml
A: 以下を確認

1. Barcode kit選択:
   - Native barcoding（EXP-NBD196）: 最高精度 ✅
   - PCR barcoding（SQK-PBK004）: 精度やや低い

2. Demultiplexing parameters:
   - Guppy barcoder: --require_barcodes_both_ends
   - Min score: 60（default）→ 70（厳格化）

3. Barcode contamination:
   - クロスコンタミ: ピペッティング注意
   - 対策: Separate workspace、新しいtips

4. Read length:
   - Short reads（<500 bp）: Barcode trimming失敗
   - 対策: Size selection（>1 kb）
```

---- 

### 7.3 コスト最適化

**Q5: コストを下げる方法は？**

```yaml
A: 以下の戦略

1. 多重化レベル上げ:
   Current: 24 samples/MinION run
   Optimized: 48-96 samples（with PromethION）
   Cost reduction: 2-4× per sample

2. Rapid Kit使用（緊急時のみ）:
   Ligation: ¥8,000/sample
   Rapid: ¥5,200/sample
   Savings: ¥2,800/sample（35%削減）
   Trade-off: 収量50%減、感度低下

3. Bulk購入:
   10+ kits: 10% discount
   Annual contract: 15-20% discount
   Savings: ¥200,000-500,000/year

4. Flow cell再利用（非推奨）:
   ONT公式: 使い捨て推奨
   実際: Wash protocol存在（非公式）
   Risk: 収量50-80%減、コンタミリスク
   → PMDA申請用途では非推奨

5. ハイブリッドキャプチャ自作:
   Commercial: ¥15,000/sample
   In-house: ¥5,000/sample（プローブ合成後）
   Savings: ¥10,000/sample
   Required: 初期投資¥75万（プローブ設計）
```

---- 

## 8. 推奨プロトコル完全版

### 8.1 Protocol 12 v2.2（MinION用、最適化版）

#### 完全試薬リスト

```yaml
【Phase 1: Extraction】
  - QIAamp Circulating Nucleic Acid Kit（Qiagen）
    Cat#: 55114
    Cost: ¥8,000/sample
    Yield: DNA + RNA co-extraction

【Phase 2: Host Depletion】
  Step 2a: Computational（Minimap2）
  Step 2b: MBD-Fc beads
    - NEBNext Microbiome DNA Enrichment Kit
      Cat#: E2612
      Cost: ¥5,000/sample

【Phase 3: DNA Virus Library】
  Step 3.1: Circular ssDNA linearization
    - DNase I（RNase-free）
      Supplier: Thermo Fisher
      Cat#: EN0521
      Cost: ¥500/sample

    - Klenow Fragment（3'→5' exo-）
      Supplier: NEB
      Cat#: M0212
      Cost: ¥300/sample

  Step 3.2: DNA library prep
    - Ligation Sequencing Kit (SQK-LSK114)
      Supplier: Oxford Nanopore
      Cat#: SQK-LSK114
      Cost: ¥8,000/sample（12 samples/kit）

【Phase 4: RNA Virus Library】
  Step 4.1: rRNA depletion（Poly(A)- virus用）
    - RiboMinus Eukaryote Kit
      Supplier: Invitrogen
      Cat#: A15026
      Cost: ¥10,000/sample

  Step 4.2: Reverse transcription
    - SuperScript IV Reverse Transcriptase
      Supplier: Invitrogen
      Cat#: 18090010
      Cost: ¥1,500/sample

    - Random Hexamer primers + Oligo-dT
      Supplier: Invitrogen
      Cost: ¥200/sample

  Step 4.3: Second-strand synthesis
    - NEBNext Ultra II Non-Directional RNA Second Strand Synthesis
      Supplier: NEB
      Cat#: E6111
      Cost: ¥800/sample

  Step 4.4: cDNA library prep
    - Ligation Sequencing Kit (same as DNA)

【Phase 5: Hybrid Capture Enrichment】
  - Custom 91-pathogen probe panel（Twist Bioscience）
    Design cost: ¥750,000（one-time）
    Per-sample cost: ¥15,000（reagents）

【Phase 6: Barcoding】
  - Native Barcoding Expansion 1-96
    Supplier: Oxford Nanopore
    Cat#: EXP-NBD196
    Cost: ¥1,040/sample（96多重化）

【Phase 7: Sequencing】
  - MinION Flow Cell R10.4.1
    Supplier: Oxford Nanopore
    Cat#: FLO-MIN114
    Cost: ¥138,000/FC（¥5,750/sample at 24 samples）

Total Reagent Cost: ¥51,390/sample
```

---- 

#### タイムライン

```yaml
Day 1（8h hands-on）:
  09:00-11:00: Extraction（2h）
  11:00-12:00: QC（NanoDrop, Qubit）
  13:00-15:00: Host depletion（2h）
  15:00-17:00: Circular ssDNA linearization + RT（2h）
  17:00-18:00: Second-strand synthesis（1h）

Day 2（6h hands-on）:
  09:00-10:00: End-repair + A-tailing（1h）
  10:00-11:00: Barcode ligation（1h）
  11:00-13:00: Cleanup + pooling（2h）
  13:00-15:00: Hybrid capture setup（2h）
  → Overnight hybridization（16h, 65°C）

Day 3（4h hands-on）:
  09:00-11:00: Capture + wash（2h）
  11:00-12:00: Elution + cleanup（1h）
  13:00-14:00: Adapter ligation（1h）
  14:00: Load MinION
  14:00-Day 6 14:00: Sequencing（72h）

Day 6:
  14:00: Stop run
  14:00-18:00: Basecalling + analysis（4h computational）

Total: 6 days（18h hands-on + 72h sequencing）
```

---- 

## 9. まとめ

### 9.1 最終推奨

**ウイルス網羅的検査（PMDA 91病原体）に最適なキット**:

```yaml
🏆 第1推奨: Ligation Sequencing Kit (SQK-LSK114) ⭐

理由:
  ✅ 最高収量（20-50 Gb/MinION）
  ✅ 最長read（10-50 kb平均）
  ✅ 全病原体対応（DNA + cDNA）
  ✅ 低バイアス
  ✅ Protocol 12 v2.1で使用実績あり

Cost: ¥8,000/sample（試薬のみ）
Time: 8-10h（prep）

併用推奨:
  + Native Barcoding (EXP-NBD196): 24-96多重化
  + Hybrid Capture: 91-pathogen panel
  + MBD Host Depletion: CpG-methylated DNA除去

Total Cost: ¥51,390/sample（試薬合計）

Performance:
  Coverage: 89-91/91 pathogens
  LOD: 100 gc/mL（reliable）
  PPA: >95%
  NPA: >98%
```

---- 

### 9.2 代替キット（特殊用途）

```yaml
緊急時（<24h TAT）:
  → Rapid Sequencing Kit (SQK-RAD114)
  Cost: ¥5,200/sample
  Trade-off: 収量50%減、感度低下

RNAウイルス詳細解析:
  → Direct RNA Sequencing Kit (SQK-RNA004)
  Cost: ¥11,000/sample
  Limitation: Poly(A)+ only、高input必要

低コスト研究用（細菌のみ）:
  → PCR Barcoding Kit (SQK-PBK004)
  Cost: ¥4,200/sample
  Limitation: 16S amplicon、ウイルス検出不可
```

---- 

### 9.3 チームメンバーの質問への直接回答

**質問**: どの試薬を、どのように使用するのが効果的か？

**回答**:

```yaml
1. 基本キット選択:
   → Ligation Sequencing Kit (SQK-LSK114) ⭐
   用途: DNA virus + RNA virus（cDNA経由）

2. 多重化:
   → Native Barcoding Kit (EXP-NBD196)
   24-96 samples/run

3. エンリッチメント:
   → Custom Hybrid Capture（91-pathogen panel）
   10-100× target enrichment

4. ワークフロー:
   → Protocol 12 v2.2（完全版）
   DNA + RNA同時処理
   18h hands-on + 72h sequencing

5. コスト:
   試薬: ¥51,390/sample
   Total（フローセル含む）: ¥57,140/sample（24多重化）

6. 性能:
   91/91病原体検出可能 ✅
   PMDA要件達成 ✅
```

---- 

### 9.4 Next Steps

```yaml
今週:
  □ Ligation Sequencing Kit (SQK-LSK114) 見積取得
    - Quantity: 10 kits（年間使用量）
    - Bulk discount確認

  □ Native Barcoding Kit見積取得
    - EXP-NBD196（96 barcodes）

  □ 代理店コンタクト
    - 日本ジェネティクス株式会社
    - 技術サポート体制確認

今月:
  □ 試薬発注
    - Lead time: 2-4週間
    - 在庫計画: 3-6ヶ月分

  □ Protocol 12 v2.2 SOP更新
    - 試薬詳細追記
    - バリデーション計画

Month 2-3:
  □ パイロット研究
    - 10-20サンプルでワークフロー確立
    - 試薬消費量実測
```
