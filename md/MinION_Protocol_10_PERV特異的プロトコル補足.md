# MinION用メタゲノム解析プロトコル
# 第10章: PERV特異的プロトコル補足

## 目次

1. [PERVの重要性と異種移植リスク](#1-pervの重要性と異種移植リスク)
2. [PERV分子生物学的特徴](#2-perv分子生物学的特徴)
3. [PERV検出戦略](#3-perv検出戦略)
4. [PERV特異的MinIONプロトコル](#4-perv特異的minionプロトコル)
5. [PERV定量と感染性評価](#5-perv定量と感染性評価)
6. [PERVサブタイプ分類](#6-pervサブタイプ分類)
7. [PERV不活化ブタの検証](#7-perv不活化ブタの検証)
8. [ヒト細胞感染性試験](#8-ヒト細胞感染性試験)
9. [Long-term モニタリング戦略](#9-long-term-モニタリング戦略)
10. [規制対応とドキュメンテーション](#10-規制対応とドキュメンテーション)

---

## 1. PERVの重要性と異種移植リスク

### 1.1 PERVとは

**豚内因性レトロウイルスの基礎**

```yaml
PERV (Porcine Endogenous Retrovirus):

  分類:
    科: Retroviridae
    属: Gammaretrovirus
    特徴: 豚ゲノムに内在化したレトロウイルス

  ゲノム構造:
    サイズ: 約9 kb
    構成: 5'-LTR-gag-pol-env-LTR-3'
    コピー数: 豚ゲノムあたり50-100コピー

  特徴:
    - 全ての豚が保有 (内因性)
    - ゲノムDNAに組み込み
    - 垂直伝播 (親から子へ)
    - 一部は複製能保持 (replication-competent)

異種移植における危険性:

  1. ヒト細胞への感染能:
     - PERV-A: ヒト細胞に感染可能
     - PERV-B: ヒト細胞に感染可能
     - PERV-C: 豚細胞特異的 (ヒトに非感染性)

  2. 組換えリスク:
     - PERV-A/C組換え体: 高感染性
     - 異種移植後のレシピエント体内で組換え可能性

  3. 長期潜伏:
     - プロウイルスとしてヒトゲノム統合可能
     - 数年〜数十年後の発症リスク

  4. 予測不能性:
     - 異種環境での挙動不明
     - 長期影響不明
```

### 1.2 PMDAガイドラインとPERV

**日本の規制要求**

```yaml
PMDAガイドライン要求事項:

  1. PERV検出の必須性:
     「異種移植用ドナー動物は、既知および未知の感染性病原体について
      検査すること。特に内因性レトロウイルスについては重点的に評価」

  2. 検出方法:
     - DNA level: プロウイルス検出
     - RNA level: 転写活性評価
     - タンパク質 level: 発現確認
     - 感染性試験: ヒト細胞への感染能評価

  3. 定量要求:
     - コピー数定量 (qPCR)
     - ウイルス粒子定量 (RT-qPCR, TCID50)

  4. サブタイプ同定:
     - PERV-A, B, C の区別
     - 組換え体検出

  5. 長期モニタリング:
     - ドナー動物の継続的監視
     - レシピエント追跡調査

PERV不活化の重要性:

  CRISPR/Cas9技術により:
    - 全PERVコピー不活化可能 (例: 69編集)
    - 本プロトコルでの検証対象:
      * 全59 PERV lociの不活化確認
      * Residual PERV activity検出
      * Off-target編集確認
```

---

## 2. PERV分子生物学的特徴

### 2.1 ゲノム構造詳細

**MinION long-read sequencingに最適な標的**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PERVゲノム構造 (約9 kb)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

5'-LTR (634 bp)
  ├─ U3 region (promoter)
  ├─ R region
  └─ U5 region

gag gene (1,500 bp)
  ├─ MA (matrix protein)
  ├─ CA (capsid protein)
  └─ NC (nucleocapsid protein)

pol gene (3,200 bp)
  ├─ PR (protease)
  ├─ RT (reverse transcriptase) ← 重要検出標的
  └─ IN (integrase)

env gene (2,000 bp)
  ├─ SU (surface protein) ← サブタイプ決定領域
  └─ TM (transmembrane protein)

3'-LTR (634 bp)
  同じ配列 (5'-LTRと)

PBS (Primer Binding Site)
PPT (Polypurine Tract)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MinION検出戦略:

  Full-length PERV sequencing:
    - 目標: 完全長9 kbリード取得
    - 利点: 全遺伝子・LTRを単一リードで解析
    - サブタイプ、組換え、変異を包括的評価

  重要領域:
    1. env gene (2 kb):
       - サブタイプ決定 (A, B, C)
       - 受容体結合領域
       - 中和エピトープ

    2. pol gene RT region (1 kb):
       - 保存領域 (universal primer設計)
       - 系統解析
       - CRISPR標的領域確認

    3. LTR (0.6 kb):
       - プロモーター活性
       - 転写レベル評価
       - 組込み部位同定

MinIONの優位性:
  - 完全長配列: サブタイプ・組換え正確判定
  - Long read: LTR-LTR junction読める
  - Direct RNA-seq: 転写活性直接評価 (オプション)
```

### 2.2 PERV多様性とコピー数

**個体・品種間変動**

```yaml
豚ゲノム内PERVコピー数:

  野生型豚:
    総コピー数: 50-100コピー/genome
    内訳:
      - Intact (完全長): 10-20コピー
      - Defective (欠損): 30-80コピー
      - Solo LTR: 多数

  Yucatan miniature pig:
    総コピー数: 約70コピー
    本プロジェクト標的: 59 PERV loci不活化

  コピー間多様性:
    配列同一性: 85-99%
    完全同一コピー: 稀
    → 各コピー個別検出必要

Replication-competent PERV:

  定義: 完全長ORF保持、感染性粒子産生可能

  頻度:
    野生型: 5-15コピー/genome
    不活化ブタ: 0コピー (目標)

  検出重要性:
    - 異種移植リスク直結
    - 不活化検証の核心
```

---

## 3. PERV検出戦略

### 3.1 多層検出アプローチ

**DNA/RNA/タンパク質レベル統合評価**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Level 1: DNA検出 (プロウイルス)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

目的: ゲノム統合PERVの検出

方法:
  1. MinION whole genome sequencing:
     - 全PERVコピー網羅的検出
     - 統合部位同定
     - 不活化部位確認

  2. Targeted sequencing:
     - PERV特異的enrichment
     - pol gene完全長シーケンス
     - env gene完全長シーケンス

  3. qPCR (定量):
     - 総PERVコピー数
     - Intact vs Defective比

判定基準:
  野生型: 50-100コピー検出
  不活化ブタ: 全標的lociで不活化配列確認
              Intact PERVゼロ

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Level 2: RNA検出 (転写活性)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

目的: PERVの転写活性評価

方法:
  1. MinION cDNA-Seq (本プロトコル):
     - 血漿cfRNA由来cDNA
     - Full-length PERV transcriptシーケンス
     - Spliced vs Unspliced RNA区別

  2. RT-qPCR:
     - gag, pol, env mRNA定量
     - サブタイプ別定量

  3. Northern blot (オプション):
     - Genomic vs Spliced RNA区別

判定基準:
  野生型: PERV RNA検出 (103-106 copies/mL血漿)
  不活化ブタ: PERV RNA検出限界以下
              または >1,000倍減少

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Level 3: タンパク質検出
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

目的: ウイルス粒子産生確認

方法:
  1. ELISA:
     - p27 Gag protein (capsid)
     - gp70 Env protein (envelope)

  2. Western blot:
     - Gag polyprotein processing確認
     - Env glycosylation確認

  3. Immunofluorescence:
     - 細胞内局在
     - Budding確認

判定基準:
  野生型: p27陽性 (>100 pg/mL)
  不活化ブタ: p27陰性 (<10 pg/mL)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Level 4: 感染性試験 (最重要)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

目的: ヒト細胞への感染能評価

方法:
  1. Co-culture assay:
     - ドナー豚細胞 + ヒト細胞 (HEK293, HT-1080)
     - 4週間培養
     - ヒト細胞からのPERV DNA検出

  2. Marker rescue assay:
     - 血漿→ヒト細胞感染
     - Nested PCR検出

  3. TCID50測定:
     - ウイルス力価定量

判定基準:
  野生型: ヒト細胞感染陽性
  不活化ブタ: ヒト細胞感染陰性
              4回反復で全て陰性
```

---

## 4. PERV特異的MinIONプロトコル

### 4.1 PERV-targeted library調製

**Enrichment戦略**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
方法1: PCR-based enrichment (推奨)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

PERV full-length amplification:

  プライマー設計:
    Forward (5'-LTR):
      5'-GCCACCCTATACTAGGAACTGAGC-3'
      位置: 5'-LTR U3 region

    Reverse (3'-LTR):
      5'-CTGACCTGCCCAGTGCTGGGATTAAAG-3'
      位置: 3'-LTR U5 region

    期待産物: 8.5-9.2 kb (full-length PERV)

  PCR条件:
    Polymerase: LongAmp Taq (NEB)
                または PrimeSTAR GXL (Takara)
    Template: cfDNA 50 ng

    反応系 (50 μL):
      Template DNA: 50 ng
      Forward primer (10 μM): 2.5 μL
      Reverse primer (10 μM): 2.5 μL
      2× LongAmp Master Mix: 25 μL
      H2O: to 50 μL

    サーマルサイクル:
      94°C - 1 min
      ↓
      30 cycles:
        94°C - 15 sec
        60°C - 30 sec
        65°C - 9 min (1 kb/min)
      ↓
      65°C - 10 min (final extension)
      4°C - hold

  精製:
    AMPure XP beads (0.6× ratio)
    溶出: 20 μL EB

  収量期待:
    野生型: 200-500 ng
    不活化ブタ: <10 ng (不活化によりPCR失敗期待)

  MinION library調製:
    SQK-LSK114使用
    通常プロトコル (第5章参照)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
方法2: CRISPR-Cas9 based enrichment
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

原理: Cas9でPERV flanking領域切断、PERV領域濃縮

  使用キット:
    - Oxford Nanopore Cas9 Sequencing Kit

  gRNA設計:
    PERV両端に2本gRNA設計
    切断→PERV領域のみ回収

  利点:
    - PCR biasなし
    - Intact/Defective両方検出
    - 統合部位情報保持

  欠点:
    - コスト高
    - 手順複雑

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
方法3: Hybrid capture enrichment
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

原理: ビオチン標識プローブでPERV配列捕捉

  プローブ設計:
    PERV全長配列カバーする120 merプローブ
    タイル状配置 (60 bp overlap)

  プロトコル:
    1. MinION library調製 (SQK-LSK114)
    2. Hybridization (65°C, 24 hours)
    3. Streptavidin beads capture
    4. 洗浄、溶出
    5. シーケンス

  利点:
    - 複数PERVコピー同時捕捉
    - Quantitative
    - 統合部位情報保持

  推奨: 研究グレード解析
```

### 4.2 PERV RNA検出プロトコル

**転写活性評価**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
cfRNA抽出 (PERV特異的考慮)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  通常プロトコル (第3章) に追加:

  RNase inhibitor強化:
    - RNaseOUT濃度2倍
    - 理由: レトロウイルスRNA不安定

  抽出直後のDNase処理:
    必須 (プロウイルスDNA混入除去)
    qPCR確認: DNA残存 <0.1%

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PERV特異的RT-PCR
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  One-step RT-PCR:
    キット: SuperScript III One-Step RT-PCR (Invitrogen)

    プライマー:
      PERV gag forward:
        5'-GTGGATCCCCGGGTCAGCTTCCGGGTAAG-3'
      PERV gag reverse:
        5'-CGGAATTCTAGAGGTCCTTTCCCCAGTTCCCT-3'
      産物: 600 bp

    反応:
      Template RNA: 100 ng
      55°C - 30 min (RT)
      94°C - 2 min
      35 cycles:
        94°C - 15 sec
        58°C - 30 sec
        68°C - 1 min
      68°C - 5 min

    判定:
      野生型: 600 bp band検出
      不活化: No band

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
cDNA全長シーケンス (MinION)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  cDNA合成:
    第6章プロトコル使用
    PERV特異的プライマー追加 (オプション):
      PERV-specific reverse primer:
        5'-CTGACCTGCCCAGTGCTGGGATTAAAG-3'
        Random primerと混合使用

  期待結果:
    野生型: Full-length PERV cDNA (9 kb) 検出
    不活化: PERV cDNA検出限界以下
```

---

## 5. PERV定量と感染性評価

### 5.1 qPCR定量プロトコル

**絶対定量とコピー数算出**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PERV DNA qPCR (総コピー数)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  標的: pol gene (保存領域)

  プライマー・プローブ:
    Forward: 5'-CTAAACCCCTCCCTGGGGTTATTAA-3'
    Reverse: 5'-GCAGAGGCCAATTGTGATGA-3'
    TaqMan probe:
      5'-FAM-ATGGCCTTGGTGGCCCCTTCT-TAMRA-3'

  Standard curve作成:
    PERV pol gene クローニング
    プラスミド希釈系列:
      10^7 - 10^1 copies/μL

  反応系 (20 μL):
    Template DNA: 5 μL (50 ng)
    Forward primer (10 μM): 0.5 μL
    Reverse primer (10 μM): 0.5 μL
    TaqMan probe (10 μM): 0.5 μL
    2× TaqMan Master Mix: 10 μL
    H2O: 3.5 μL

  サーマルサイクル:
    50°C - 2 min (UDG activation)
    95°C - 10 min
    40 cycles:
      95°C - 15 sec
      60°C - 1 min (data collection)

  コピー数計算:
    Standard curveからCt→copy数変換
    Copy/genome = (Total copies) / (Input genome数)

    Input genome数 = (DNA ng × 6.022×10^23) / (Genome size × 10^9 × 650)
    豚genome: 2.7×10^9 bp

  判定基準:
    野生型: 50-100 copies/genome
    不活化ブタ: <1 copy/genome (検出限界以下)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PERV RNA RT-qPCR (転写活性)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  標的: gag mRNA

  One-step RT-qPCR:
    キット: QuantiFast SYBR Green RT-PCR Kit

    プライマー:
      Forward: 5'-GTGGATCCCCGGGTCAGCTTCCGGGTAAG-3'
      Reverse: 5'-CGGAATTCTAGAGGTCCTTTCCCCAGTTCCCT-3'

    Standard: In vitro transcribed PERV RNA
              10^8 - 10^2 copies/μL

    反応:
      Template RNA: 5 μL (50 ng)
      Primer mix: 各0.5 μL
      2× RT-PCR Mix: 10 μL
      RT mix: 0.2 μL

    サーマルサイクル:
      50°C - 10 min (RT)
      95°C - 5 min
      40 cycles:
        95°C - 10 sec
        60°C - 30 sec (data collection)

    Melt curve: 60-95°C, 0.5°C/sec

  コピー数計算:
    copies/mL血漿 = (measured copies) × (溶出量/RNA input) × (1/血漿量)

  判定基準:
    野生型: 10^3 - 10^6 copies/mL
    不活化: <10 copies/mL (検出限界以下)
```

### 5.2 ウイルス力価測定

**感染性粒子定量**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TCID50 assay
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  原理:
    50% Tissue Culture Infectious Dose
    = 細胞の50%を感染させるウイルス量

  プロトコル:

    1. サンプル調製:
       血漿: 0.22 μm filter滅菌
       10倍段階希釈: 10^-1 - 10^-8

    2. 細胞準備:
       HEK293T cells (ヒト腎細胞)
       96-well plate: 1×10^4 cells/well
       各希釈度: 8 wells反復

    3. 感染:
       希釈ウイルス100 μL添加/well
       37°C, 5% CO2, 7日間培養

    4. 継代:
       Week 1: 培地交換、細胞継代
       Week 2-4: さらに3回継代

    5. PERV検出:
       各wellの細胞からDNA抽出
       PERV pol qPCR
       Positive wells数カウント

    6. TCID50計算:
       Reed-Muench法
       TCID50/mL = 10^(x) × dilution factor

  判定基準:
    野生型: 10^3 - 10^5 TCID50/mL
    不活化: <10 TCID50/mL (陰性)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Plaque assay (代替法)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  利点: 直接可視化

  プロトコル:
    1. HEK293T単層培養 (6-well plate)
    2. ウイルス希釈液接種 (10^-1 - 10^-6)
    3. 1時間吸着
    4. Overlay medium (0.5% agarose)
    5. 7-14日培養
    6. Crystal violet染色
    7. Plaque counting

  判定:
    Plaque Forming Units (PFU)/mL計算

  注意: PERVはplaque形成遅い・不明瞭
       → TCID50推奨
```

---

## 6. PERVサブタイプ分類

### 6.1 env gene配列解析

**サブタイプ決定プロトコル**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PERV env遺伝子増幅
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  プライマー設計 (Universal):
    env-F: 5'-TGGGCAACTTGTGCATCATCAATC-3'
    env-R: 5'-GGATGTGCATTGGTATTTCAAGTT-3'
    産物: 2.1 kb (env全長)

  PCR条件:
    Template: cfDNA 100 ng
    Polymerase: PrimeSTAR GXL

    反応系 (50 μL):
      Template: 100 ng
      Primer各: 0.4 μM
      2× GXL Buffer: 25 μL

    サーマルサイクル:
      98°C - 1 min
      30 cycles:
        98°C - 10 sec
        60°C - 15 sec
        68°C - 2 min
      68°C - 5 min

  精製: AMPure XP (0.8× ratio)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MinION env sequencing
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ライブラリ調製:
    SQK-LSK114使用
    Input: 200 ng PCR product

  シーケンス:
    Flowcell: 1/4使用 (multiplexing可能)
    ラン時間: 6-12時間 (2 kb短鎖のため)
    期待reads: 50,000-100,000 env reads

  Basecalling: Duplex mode (Q30+)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
サブタイプ分類基準
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  受容体結合領域 (RBD) 配列比較:

  PERV-A:
    RBD配列: WDYPL motif保存
    受容体: Hamster PaR-1 (ヒト細胞で発現)
    ヒト感染性: 陽性
    Reference: GenBank AF038599

  PERV-B:
    RBD配列: WDYQL motif
    受容体: Hamster PaR-2 (ヒト細胞で発現)
    ヒト感染性: 陽性
    Reference: GenBank AF038600

  PERV-C:
    RBD配列: WDNKP motif
    受容体: 豚特異的受容体
    ヒト感染性: 陰性
    Reference: GenBank AF038601

  分類手順:
    1. MinION read alignment to references
    2. RBD motif抽出 (position 120-125)
    3. Motif pattern matching
    4. サブタイプ割り当て

  組換え体検出:
    - PERV-A/C recombinant: 高ヒト感染性
    - Mosaic pattern検出
    - Breakpoint同定

  バイオインフォマティクス:
    # Minimap2 alignment
    minimap2 -ax map-ont \
      perv_env_refs.fa \
      env_reads.fastq | \
      samtools view -bS - > env_aligned.bam

    # Subtype calling script
    python perv_subtype_classifier.py \
      --bam env_aligned.bam \
      --out subtype_report.txt

  出力例:
    Sample: Pig001_plasma
    Total PERV env reads: 8,523

    Subtype distribution:
      PERV-A: 3,245 (38.1%)
      PERV-B: 2,876 (33.7%)
      PERV-C: 2,155 (25.3%)
      Recombinant A/C: 247 (2.9%)

    Risk assessment:
      Human-tropic subtypes (A+B+A/C): 74.7%
      High-risk recombinants: 2.9%
```

---

## 7. PERV不活化ブタの検証

### 7.1 CRISPR編集部位確認

**全59 loci検証プロトコル**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Whole genome sequencing approach
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  目的: 全PERVコピーの不活化確認

  MinION ultra-long read sequencing:
    Input: High molecular weight genomic DNA
           (cfDNAでなく、組織DNA使用)

    抽出: Nanobind CBB Big DNA Kit
    Library: SQK-LSK114
    Flowcells: 2-3枚
    目標カバレッジ: 30× whole genome

  PERV loci mapping:
    1. Alignment to Sus scrofa reference + PERV loci annotation
    2. 各59 loci位置の読み取り
    3. CRISPR編集痕跡確認:
       - Indel (insertion/deletion)
       - 完全欠失
       - Substitution

  判定:
    合格: 全59 lociで編集確認
          Intact PERVゼロ
    不合格: 1つでもIntact PERVあり

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Targeted amplicon sequencing (検証用)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  各CRISPR標的部位のPCR増幅:

  プライマー設計:
    各locus flanking領域にプライマー設計
    59セット準備

  Multiplex PCR:
    全59アンプリコン同時増幅
    Pooling後MinIONシーケンス

  編集効率算出:
    Wild-type reads: 編集なし
    Edited reads: indel/deletion

    Editing efficiency = (Edited / Total) × 100%

  目標: 全lociで >99% editing efficiency

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Off-target編集確認
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  予測off-target部位:
    - Cas-OFFinder使用
    - ミスマッチ<3 bpの配列リストアップ

  確認:
    Whole genome sequencingデータから
    予測off-target部位の配列抽出

    判定:
      合格: 編集痕跡なし
      要調査: 編集検出 (影響評価必要)
```

### 7.2 Residual PERV活性スクリーニング

**微量残存活性検出**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
超高感度PCR
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Nested PCR (2段階増幅):

    1st PCR:
      Template: cfDNA 500 ng (通常の10倍)
      Primer: PERV pol outer primers
      Cycles: 30

    2nd PCR:
      Template: 1st PCR product 1 μL
      Primer: PERV pol inner primers (nested)
      Cycles: 30

    検出感度: <1 copy/500 ng DNA
              = <0.001 copies/genome

  判定:
    野生型: Strong band
    不活化: No band (60反復で全て陰性確認)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Digital droplet PCR (ddPCR)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  絶対定量・超高感度:

  機器: QX200 Droplet Digital PCR System (Bio-Rad)

  プロトコル:
    1. Droplet生成: 20,000 droplets/sample
    2. PCR (PERV pol TaqMan assay)
    3. Droplet読み取り
    4. Poisson統計で絶対濃度算出

  検出限界: 0.0001 copies/genome

  判定:
    不活化確認: <0.001 copies/genome
                = 実質ゼロ

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Long-term monitoring
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  頻度: 月1回 × 12ヶ月

  サンプル:
    - 血漿 (cfDNA/RNA)
    - 唾液
    - 尿

  測定:
    - PERV DNA qPCR
    - PERV RNA RT-qPCR
    - Nested PCR
    - Co-culture assay (3ヶ月ごと)

  判定:
    全timepoint陰性: 不活化成功
    1回でも陽性: 詳細調査、動物除外
```

---

## 8. ヒト細胞感染性試験

### 8.1 Co-culture assay

**Gold standardテスト**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
プロトコル
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  細胞準備:

    Donor cells (ブタ由来):
      - PK-15 (豚腎細胞)
      - または primary pig cells
      - PERV産生能確認済み

    Target cells (ヒト由来):
      - HEK293T (ヒト腎細胞)
      - HT-1080 (ヒト線維肉腫)
      - 両方使用推奨

  Co-culture設定:

    1. Target cells播種:
       6-well plate
       5×10^5 cells/well
       前日準備

    2. Donor cells添加:
       5×10^5 cells/well
       1:1 ratio

    3. 共培養:
       37°C, 5% CO2
       4週間

    4. 継代:
       週1回
       Target cells選択的継代:
         - Trypsinization (0.05%)
         - 豚細胞は剥がれにくい
         - ヒト細胞のみ回収

    5. Week 4でPERV検出:
       Target cells DNA抽出
       PERV pol qPCR
       Nested PCR

  Positive control:
    PERV感染HEK293T細胞 (既知陽性)

  Negative control:
    Target cellsのみ (ブタ細胞なし)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
判定基準
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  陽性 (感染あり):
    qPCR Ct <35
    かつ Nested PCR band陽性
    → ヒト細胞感染性あり (異種移植不可)

  陰性 (感染なし):
    qPCR Ct >38 (検出限界以下)
    かつ Nested PCR band陰性
    → ヒト細胞感染性なし

  判定信頼性向上:
    - 4反復実施
    - 全て陰性で初めて「陰性」判定
    - 1つでも陽性: 動物除外
```

### 8.2 Pseudotype virus assay

**env機能評価**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
原理
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  PERVのenvタンパク質を持つが
  ゲノムはreporter gene (GFP, luciferase)

  → 1回のみの感染検出
  → 安全性高い

プロトコル:

  1. Pseudovirus作製:
     HEK293T transfection:
       - PERV env発現プラスミド
       - MLV gag-pol プラスミド
       - Luciferase reporter プラスミド

     48時間後: 上清回収 (pseudovirus含有)

  2. Target cells感染:
     HEK293T, HT-1080播種
     Pseudovirus添加

  3. 48時間後: Luciferase測定

  判定:
    RLU (Relative Light Units) >2× background:
      → PERV env機能的、ヒト細胞侵入可能

    RLU ≤ background:
      → PERV env非機能的または不在

  サブタイプ別評価:
    PERV-A env → RLU高い
    PERV-B env → RLU中程度
    PERV-C env → RLU低い (ヒト細胞非感染性)
```

---

## 9. Long-term モニタリング戦略

### 9.1 ドナー動物モニタリング

**継続的PERV評価**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
モニタリングスケジュール
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Pre-breeding検査 (繁殖前):
    - Whole genome sequencing
    - 全PERV loci編集確認
    - Co-culture assay (4週間)
    - 判定: 全て陰性で繁殖許可

  Monthly検査 (生涯継続):
    - 血漿cfDNA/RNA抽出
    - PERV qPCR/RT-qPCR
    - 唾液・尿PERV PCR

  Quarterly検査 (3ヶ月ごと):
    - Co-culture assay
    - 血清p27 ELISA

  Annual検査 (年1回):
    - Whole genome re-sequencing
    - PERV loci再確認
    - Pseudotype virus assay

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
データベース管理
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  個体別追跡:
    Animal ID: Unique identifier
    全検査結果の時系列記録

    Database schema:
      - Animal_ID
      - Sampling_Date
      - Sample_Type (blood, saliva, urine)
      - PERV_DNA_qPCR (copies/genome)
      - PERV_RNA_qPCR (copies/mL)
      - Co_culture (Positive/Negative)
      - Notes

  Alert system:
    PERV検出時:
      → 即座通知
      → 隔離
      → 詳細調査
      → 繁殖停止

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
世代追跡
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  F0 (CRISPR編集世代):
    最厳格モニタリング
    Off-target編集確認

  F1, F2... (子孫):
    遺伝的安定性確認
    PERV loci遺伝パターン
    新規PERV出現監視

  目標:
    F3世代以降も全PERV陰性維持
    10年以上のデータ蓄積
```

### 9.2 レシピエント追跡 (将来の臨床応用)

**ヒト患者モニタリング計画**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
移植前検査
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ベースライン確立:
    - ヒト血液PERV PCR (陰性確認)
    - 抗PERV抗体 (陰性確認)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
移植後モニタリング
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  頻度:
    Week 1-4: 週1回
    Month 2-12: 月1回
    Year 2-5: 3ヶ月ごと
    Year 6-: 6ヶ月ごと (生涯)

  検査項目:
    1. PERV DNA PCR (血液PBMC)
    2. PERV RNA RT-PCR (血漿)
    3. 抗PERV抗体 (ELISA)
    4. PBMC-pig cell co-culture (年1回)

  陽性時対応:
    - 確認検査 (重複測定)
    - 詳細配列解析 (豚由来確認)
    - 感染源特定
    - 臨床経過観察
    - Public health報告

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
家族・接触者調査
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  対象:
    - 配偶者
    - 同居家族
    - 医療従事者

  検査:
    年1回PERV PCR

  水平伝播監視:
    PERV陽性者出現 → 詳細疫学調査
```

---

## 10. 規制対応とドキュメンテーション

### 10.1 PMDA提出用データパッケージ

**包括的PERV評価レポート**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
必須提出データ
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. Donor animal characterization:
     - 個体情報 (ID, 系統, 年齢, 性別)
     - CRISPR編集詳細 (標的配列, guide RNA)
     - Whole genome sequencing data
     - PERV loci annotation

  2. PERV genomic analysis:
     - 全PERV loci配列 (FASTA format)
     - Editing efficiency (各locus)
     - Off-target analysis結果
     - Copy number (qPCR)

  3. PERV expression analysis:
     - RNA-seq data (FastQ, BAM)
     - RT-qPCR結果 (生データ + 解析)
     - Northern blot (実施時)

  4. PERV protein analysis:
     - p27 ELISA結果 (時系列)
     - Western blot images
     - Immunofluorescence images

  5. Infectivity data:
     - Co-culture assay結果 (4反復 × 3 timepoints)
     - Pseudotype virus assay
     - TCID50測定 (実施時)

  6. Subtyping data:
     - env配列解析
     - Subtype distribution
     - Recombinant検出結果

  7. Longitudinal data:
     - Monthly qPCR results (12ヶ月)
     - Quarterly co-culture (4回)
     - Trend analysis graphs

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
レポート構成
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  PERV Risk Assessment Report for Xenotransplantation

  Executive Summary (2-3 pages)
    - PERV不活化の成功確認
    - ヒト細胞感染性陰性
    - リスク評価結論

  1. Introduction
     1.1 Background and rationale
     1.2 Regulatory framework (PMDA guideline)
     1.3 Study objectives

  2. Materials and Methods
     2.1 Donor animal description
     2.2 CRISPR/Cas9 editing strategy
     2.3 Genomic analysis methods
         - Whole genome sequencing
         - MinION long-read sequencing
     2.4 Expression analysis methods
         - RNA extraction (cfRNA)
         - cDNA synthesis
         - RT-qPCR
     2.5 Infectivity assay methods
         - Co-culture protocol
         - Pseudotype virus assay
     2.6 Quality control and validation

  3. Results
     3.1 PERV genomic characterization
         - Table: All 59 PERV loci status
         - Figure: Editing efficiency
         - Figure: Genome browser view

     3.2 PERV expression analysis
         - Table: qPCR Ct values (time-series)
         - Figure: RT-qPCR results
         - Conclusion: No detectable PERV RNA

     3.3 PERV protein analysis
         - Table: p27 ELISA results
         - Figure: Western blot
         - Conclusion: No PERV protein

     3.4 PERV infectivity analysis
         - Table: Co-culture results (4 replicates)
         - Figure: Target cell PERV PCR
         - Conclusion: No human cell infection

     3.5 PERV subtyping
         - Table: env sequence analysis
         - Figure: Phylogenetic tree
         - Conclusion: All subtypes inactivated

     3.6 Longitudinal monitoring
         - Figure: 12-month trend
         - Table: All monthly results
         - Conclusion: Stable PERV-negative status

  4. Risk Assessment
     4.1 Probability of PERV transmission
         - Genomic: Negligible (all loci inactivated)
         - Expression: Negligible (no RNA/protein)
         - Infectivity: Negligible (all assays negative)
         Overall: <0.001% probability

     4.2 Severity if transmission occurs
         - Based on literature review
         - Worst-case scenarios

     4.3 Overall risk level
         - Risk matrix
         - Conclusion: Acceptable for clinical trial

  5. Discussion
     5.1 Comparison with wild-type pigs
     5.2 Limitations of the study
     5.3 Future monitoring plan

  6. Conclusions
     PERV inactivation successful
     Human infection risk minimized
     Suitable for xenotransplantation

  7. References

  Appendices
     A. Raw data tables
     B. Sequence alignments
     C. Statistical analysis
     D. SOPs
     E. QC records (ALCOA+)
     F. Reagent lot numbers
     G. Equipment calibration records

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
データ管理とトレーサビリティ
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Electronic Lab Notebook (ELN):
    全実験の詳細記録
    リアルタイム入力
    署名・承認ワークフロー

  LIMS (Laboratory Information Management System):
    サンプル追跡
    試薬ロット管理
    機器使用記録

  Data storage:
    Raw data: 本番サーバー + バックアップ
    Processed data: 解析サーバー
    Reports: Document management system
    保存期間: 最低15年 (規制要求)

  Audit trail:
    全データ変更履歴
    誰が、いつ、何を変更
    変更理由記録
```

### 10.2 ALCOA+準拠PERV検査記録

**規制監査対応**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PERV検査記録フォーマット
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  PERV Detection Record #: PERV-YYYYMMDD-001

  【Attributable (帰属性)】
    実施者: [フルネーム] _______ 署名: ______ 日時: ______
    確認者: [フルネーム] _______ 署名: ______ 日時: ______
    承認者: [フルネーム] _______ 署名: ______ 日時: ______

  【Legible (判読性)】
    全て印刷体またはタイプ入力
    修正: 二重線+署名+日時

  【Contemporaneous (同時性)】
    各ステップ実施と同時に記録
    タイムスタンプ: YYYY-MM-DD HH:MM:SS

  【Original (原本性)】
    1st entry: このシート
    Copy: 禁止 (スキャンPDFは可)

  【Accurate (正確性)】
    機器から直接転記
    計算式記載
    単位明記

  【Complete (完全性)】
    全項目記入
    N/A: 該当なし明記
    逸脱: 詳細記録

  【Consistent (一貫性)】
    フォーマット統一
    命名規則遵守
    単位統一

  【Enduring (耐久性)】
    耐久性インク使用
    専用バインダー保管
    電子backup

  【Available (利用可能性)】
    監査時即座提供可能
    索引作成
    検索容易

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
記録保管
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  物理的記録:
    場所: 温湿度管理された書庫
    期間: 15年間
    アクセス: 権限者のみ

  電子記録:
    システム: 21 CFR Part 11準拠EDMS
    バックアップ: 3箇所 (本番+バックアップ+オフサイト)
    セキュリティ: 暗号化、アクセスログ

  監査準備:
    年1回内部監査
    外部監査対応訓練
    記録の完全性定期確認
```

---

## まとめ

本章では、異種移植における最重要病原体PERVの包括的検出・評価プロトコルを詳述しました。

**重要ポイント:**

1. **PERV検出の多層アプローチ**
   - DNA (ゲノム), RNA (転写), タンパク質, 感染性
   - MinION long-readの威力: 完全長配列、サブタイプ、組換え検出

2. **PERV不活化ブタの厳格な検証**
   - 全59 PERV loci編集確認
   - Residual活性の超高感度検出
   - Off-target編集確認

3. **ヒト細胞感染性試験**
   - Co-culture assay (gold standard)
   - 4反復全て陰性で初めて合格
   - Pseudotype virusでenv機能評価

4. **長期モニタリング**
   - ドナー動物: 月次検査、生涯追跡
   - レシピエント(将来): 生涯PERV監視
   - データベース管理とAlert system

5. **規制対応**
   - PMDA提出用包括的レポート
   - 15年間のデータ保管
   - ALCOA+準拠記録

6. **MinION活用の優位性**
   - Full-length PERV sequencing (9 kb)
   - サブタイプ・組換え検出
   - 統合部位同定
   - 編集部位確認

**全プロトコル完成**

第1章〜第10章により、MinION用メタゲノム解析の完全なプロトコル体系が確立されました。PMDA 91病原体検出、特にPERV評価に対応した、規制準拠の包括的システムです。

---

**文書情報**
- 作成日: 2025-03-08
- バージョン: 1.0
- 作成者: MinIONメタゲノム解析プロトコル開発チーム
- 承認者: [承認者名]
- 次回改訂予定: 2025-09-08
