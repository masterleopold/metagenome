# PERV検証方法：詳細プロトコル

**対象**: ゲノム編集PERVフリーブタの包括的検証
**規制当局**: FDA/PMDA準拠

---- 

## レベル1：DNA/遺伝子レベルの検証

### 1.1 全ゲノムシーケンス（WGS）

#### 目的
CRISPR/Cas9によるPERV遺伝子座の編集成功を確認し、オフターゲット変異を検出する。

#### プロトコル

**サンプル準備**:
```
1. 血液または組織サンプル採取（10-20 ml血液または1-2 g組織）
2. ゲノムDNA抽出（Qiagen DNeasy Blood & Tissue Kit）
3. DNA品質確認:
   - 濃度: >50 ng/μl
   - 純度: A260/280 = 1.8-2.0
   - 完全性: ゲル電気泳動で高分子量バンド確認
```

**ライブラリー調製と配列決定**:
```
1. ライブラリー調製: Illumina TruSeq DNA PCR-Free
2. カバレッジ: 30-50×（全ゲノム）
3. リード長: 150 bp paired-end
4. シーケンサー: Illumina NovaSeq 6000またはNextSeq 2000
```

**バイオインフォマティクス解析**:
```python
# 解析パイプライン
pipeline = {
    "alignment": "BWA-MEM (Sscrofa11.1参照ゲノム)",
    "variant_calling": "GATK HaplotypeCaller",
    "PERV_detection": {
        "method": "BLAST検索 + de novo assembly",
        "targets": "全62 PERV遺伝子座",
        "threshold": "99% identity, >500 bp alignment"
    },
    "off_target_analysis": {
        "tool": "Cas-OFFinder + CRISPOR",
        "predicted_sites": "上位100オフターゲット候補",
        "validation": "各候補サイトのSanger sequencing確認"
    }
}
```

**合格基準**:
- ✅ 全62 PERV遺伝子座でCRISPR編集確認（インデル、欠失）
- ✅ 野生型PERV配列の不在
- ✅ 予測オフターゲット部位に変異なし（または影響のない領域）
- ✅ 機能遺伝子の破壊的変異なし

**タイムライン**: 3-5日（ライブラリー調製1日、シーケンス2-3日、解析1日）
**コスト**: ¥180,000/サンプル
**設備**: 次世代シーケンサー（大学コアファシリティ利用推奨）

---- 

### 1.2 Droplet Digital PCR（ddPCR）

#### 目的
PERVプロウイルスのコピー数を絶対定量する。

#### プロトコル

**プライマー/プローブ設計**:
```
Target: PERV pol遺伝子（PERV-A/B/C共通領域）
Forward primer: 5'-CTGACGTTAGTGATGCTGCC-3'
Reverse primer: 5'-AGCGAAGAGTAAGGGATAGC-3'
Probe (FAM): 5'-/56-FAM/TCCCACTGGCTCGAGGGATCT/3BHQ_1/-3'

Reference gene: Sus scrofa β-actin
Forward primer: 5'-TGTTGCCCTAGACTTCGAGC-3'
Reverse primer: 5'-GGACTTCGAGCAGGAGATGG-3'
Probe (HEX): 5'-/5HEX/ATCCACACAGAGTACTTGCGCTCA/3IABkFQ/-3'
```

**ddPCR反応**:
```
1. 反応液調製（20 μl）:
   - ddPCR Supermix for Probes: 10 μl
   - Forward primer (18 μM): 1 μl
   - Reverse primer (18 μM): 1 μl
   - Probe (5 μM): 1 μl
   - Template DNA (10 ng/μl): 2 μl
   - Nuclease-free water: 5 μl

2. ドロップレット生成:
   - QX200 Droplet Generator使用
   - 70 μl droplet generation oil添加
   - 約20,000ドロップレット/ウェル生成

3. PCR増幅:
   - 95°C 10分（初期変性）
   - [94°C 30秒 → 60°C 1分] × 40サイクル
   - 98°C 10分（最終伸長）
   - 4°C 保持
   - Ramp rate: 2°C/秒

4. ドロップレット読み取り:
   - QX200 Droplet Readerで測定
   - FAM陽性ドロップレット数カウント
   - Poisson統計でコピー数計算
```

**データ解析**:
```python
# コピー数計算
def calculate_copy_number(positive_droplets, total_droplets, volume_per_droplet):
    """
    Poisson統計を使用したコピー数計算
    """
    import math

    # Poisson分布: λ = -ln(1 - p)
    # p = positive_droplets / total_droplets
    if positive_droplets == 0:
        return 0, 0  # コピー数0、95% CI上限

    fraction_positive = positive_droplets / total_droplets
    lambda_value = -math.log(1 - fraction_positive)

    # コピー数/μl
    copies_per_ul = lambda_value / volume_per_droplet

    # 95% 信頼区間（Poisson近似）
    ci_low = lambda_value - 1.96 * math.sqrt(lambda_value / total_droplets)
    ci_high = lambda_value + 1.96 * math.sqrt(lambda_value / total_droplets)

    return copies_per_ul, (ci_low / volume_per_droplet, ci_high / volume_per_droplet)

# 野生型ブタ: 約62コピー/ゲノム
# PERV-KOブタ: 0コピー（検出限界: <0.5コピー/ゲノム）
```

**合格基準**:
- ✅ PERVコピー数: **0コピー/ゲノム**
- ✅ β-actinコピー数: 2コピー/ゲノム（ディプロイド確認）
- ✅ 検出限界: \<0.5コピー/ゲノム（3連の平均）

**タイムライン**: 4-6時間
**コスト**: ¥25,000/サンプル
**設備**: Bio-Rad QX200 ddPCRシステム

---- 

### 1.3 標的Sangerシーケンシング

#### 目的
主要なCRISPR標的部位の詳細な配列確認。

#### プロトコル

**PCR増幅**:
```
1. 各PERV遺伝子座のフランキング領域をPCR増幅
2. プライマー設計: 編集部位から上流/下流200-300 bp
3. アンプリコンサイズ: 400-800 bp
4. PCR条件: 標準3ステップPCR（94°C変性、58-62°C結合、72°C伸長）
```

**配列決定**:
```
1. PCR産物精製: ExoSAP-ITまたはゲル精製
2. Sanger sequencing: ABI 3730xlまたは同等品
3. 両方向シーケンス（forward + reverse）
```

**解析**:
```
1. CRISPR標的部位のアライメント
2. インデル（挿入/欠失）の同定
3. フレームシフト変異の確認
4. アレル頻度の定量（ヘテロ接合性チェック）
```

**合格基準**:
- ✅ 全CRISPR標的部位でインデルまたは欠失
- ✅ ホモ接合性（両アレルで編集）
- ✅ フレームシフトによる遺伝子破壊

**タイムライン**: 2-3日
**コスト**: ¥5,000-10,000/標的部位（合計¥45,000で主要5箇所）

---- 

## レベル2：タンパク質レベルの検証

### 2.1 Western Blot（最重要）

#### 目的
PERV Gag、Pol、Envタンパク質の不在を確認する。

#### プロトコル

**サンプル準備**:
```
1. 組織サンプル:
   - 腎臓、心臓、肝臓、末梢血単核球（PBMC）
   - 各50-100 mg組織または10^7細胞

2. タンパク質抽出:
   - RIPA bufferでライセート作成
   - プロテアーゼ阻害剤カクテル添加
   - 超音波処理または機械的破砕
   - 4°C、12,000 gで10分遠心
   - 上清回収

3. タンパク質定量:
   - BCA assayまたはBradford assay
   - 濃度調整: 2-5 mg/ml
```

**SDS-PAGE電気泳動**:
```
1. サンプル調製:
   - タンパク質20-50 μg/レーン
   - 4× Laemmli buffer添加
   - 95°C、5分加熱（変性）

2. ゲル電気泳動:
   - 12% polyacrylamide gel（分離ゲル）
   - 4% stacking gel（濃縮ゲル）
   - 電気泳動条件: 80V（stacking）→ 120V（分離）、1-1.5時間
   - Protein ladder使用（10-250 kDa）
```

**転写（Transfer）**:
```
1. PVDF membrane準備:
   - メタノールで30秒活性化
   - Transfer bufferで平衡化

2. Semi-dry transfer:
   - 20V、30-45分（タンパク質サイズ依存）
   - または湿式transfer: 100V、1時間

3. 転写確認:
   - Ponceau S染色で全体転写確認
   - TBS-T（0.1% Tween-20）で洗浄
```

**免疫検出**:
```
1. ブロッキング:
   - 5% skim milkまたは5% BSA in TBS-T
   - 室温1時間または4°C一晩

2. 一次抗体反応:
   抗体リスト:
   - Anti-PERV p27 Gag (rabbit polyclonal, 1:1000)
   - Anti-PERV gp70 Env (mouse monoclonal, 1:500)
   - Anti-PERV p15E TM (rabbit polyclonal, 1:1000)
   - Loading control: Anti-β-actin (1:5000)

   反応条件: 4°C、一晩振とう

3. 洗浄:
   - TBS-T × 3回、各10分

4. 二次抗体反応:
   - HRP-conjugated anti-rabbit IgG (1:5000)
   - HRP-conjugated anti-mouse IgG (1:5000)
   - 室温1時間

5. 洗浄:
   - TBS-T × 3回、各10分

6. 化学発光検出:
   - ECL substrate添加
   - ChemiDoc imaging systemで検出
   - 露光時間: 10秒-5分（シグナル強度依存）
```

**陽性/陰性コントロール**:
```
陽性コントロール:
- PK-15細胞ライセート（豚腎臓由来、PERV産生細胞株）
- 野生型ブタ組織ライセート

陰性コントロール:
- HEK-293細胞ライセート（ヒト細胞、PERV陰性）
- 一次抗体なし（二次抗体のみ）

結果解釈:
野生型ブタ: p27 (27 kDa), gp70 (70 kDa), p15E (15 kDa)で強いバンド
PERV-KOブタ: 全てのPERVタンパク質でバンド不在
β-actin (42 kDa): 全サンプルで検出（ローディングコントロール）
```

**合格基準**:
- ✅ p27、gp70、p15Eタンパク質の**完全な不在**
- ✅ 陽性コントロールで明確なバンド
- ✅ β-actinが全サンプルで均一
- ✅ 複数組織（3種類以上）で確認

**タイムライン**: 2-3日
**コスト**: ¥50,000/動物（複数組織）

---- 

### 2.2 質量分析プロテオミクス（確認用）

#### 目的
バイアスのない方法でPERVペプチドの不在を確認する。

#### プロトコル

**サンプル準備**:
```
1. タンパク質抽出（Western Blotと同様）
2. タンパク質濃度調整: 100 μg
3. 還元・アルキル化:
   - DTT（10 mM）で37°C、30分
   - Iodoacetamide（55 mM）で室温、暗所30分
```

**トリプシン消化**:
```
1. クロロホルム/メタノール沈殿でバッファー交換
2. Trypsin添加（1:50 w/w比）
3. 37°C、一晩消化
4. Formic acid（1%）で反応停止
5. C18 ZipTipで脱塩
```

**LC-MS/MS分析**:
```
1. 液体クロマトグラフィー:
   - C18カラム（75 μm × 150 mm、3 μm粒子径）
   - 溶媒A: 0.1% formic acid in water
   - 溶媒B: 0.1% formic acid in acetonitrile
   - グラジエント: 5-35% B、90分

2. 質量分析:
   - Orbitrap Fusion Lumos または同等品
   - Positive ion mode
   - Full MS scan: 350-1800 m/z, resolution 120,000
   - MS/MS: Top 20 most intense ions, HCD fragmentation
```

**データベース検索**:
```python
# 検索パラメータ
search_params = {
    "database": "Sus scrofa UniProt + PERV proteins",
    "enzyme": "Trypsin",
    "missed_cleavages": 2,
    "peptide_tolerance": "10 ppm",
    "fragment_tolerance": "0.02 Da",
    "modifications": {
        "fixed": "Carbamidomethyl (C)",
        "variable": "Oxidation (M), Acetylation (Protein N-term)"
    },
    "FDR": "1% (PSM level)"
}

# PERV特異的ペプチド（検出対象）
perv_peptides = [
    "QNRPGPGR",  # p27 Gag
    "GIVQQQNNLLR",  # gp70 Env
    "LQARVLAVERYLK",  # p15E TM
    # ... 他20種類のPERV特異的ペプチド
]
```

**合格基準**:
- ✅ PERV特異的ペプチドの**完全な不在**
- ✅ ブタタンパク質の正常検出（5,000-10,000 proteins）
- ✅ 3連サンプルで一貫性

**タイムライン**: 1-2週間（サンプル調製2-3日、MS測定1日、解析3-5日）
**コスト**: ¥150,000-300,000/サンプル

---- 

### 2.3 免疫組織化学（IHC）

#### 目的
組織切片におけるPERVタンパク質の空間分布を可視化する。

#### プロトコル

**組織準備**:
```
1. 組織固定:
   - 10% neutral buffered formalin、24-48時間
   - または凍結切片（OCT compound包埋、-80°C保存）

2. パラフィン包埋（FFPE）:
   - 段階的エタノール脱水
   - キシレン透徹
   - パラフィン浸透・包埋

3. 薄切:
   - ミクロトームで4-5 μm切片作成
   - スライドガラスにマウント
```

**染色プロトコル（DAB法）**:
```
1. 脱パラフィン・再水和:
   - キシレン × 3回、各5分
   - 段階的エタノール（100% → 95% → 70%）
   - 蒸留水

2. 抗原賦活化:
   - Citrate buffer (pH 6.0)中でマイクロ波加熱
   - 95-100°C、15-20分
   - 室温まで冷却

3. 内因性ペルオキシダーゼブロック:
   - 3% H2O2 in methanol、10分

4. ブロッキング:
   - 10% normal goat serum in PBS、1時間

5. 一次抗体反応:
   - Anti-PERV p27、gp70、p15E（各1:100-1:500希釈）
   - 4°C、一晩

6. 洗浄:
   - PBS × 3回、各5分

7. 二次抗体反応:
   - HRP-conjugated secondary antibody (1:500)
   - 室温、30分

8. DAB発色:
   - DAB substrate添加
   - 顕微鏡下で発色モニター（1-5分）
   - 水で反応停止

9. 対比染色:
   - Hematoxylin、1分
   - 水洗、分別

10. 脱水・封入:
    - エタノール → キシレン
    - Mounting medium で封入
```

**顕微鏡観察**:
```
評価項目:
1. 染色強度:
   - 0: 染色なし
   - 1+: 弱い染色
   - 2+: 中等度染色
   - 3+: 強い染色

2. 染色パターン:
   - 細胞質染色
   - 細胞膜染色
   - 核染色（非特異的）

3. 陽性細胞割合:
   - 0%: 陰性
   - <10%: 散在性
   - 10-50%: 局所性
   - >50%: びまん性

野生型ブタ: 2-3+、びまん性（特に腎臓、リンパ組織）
PERV-KOブタ: 0、全組織で陰性
```

**合格基準**:
- ✅ 全組織切片でPERV染色陰性（スコア0）
- ✅ 陽性コントロール（野生型）で明確な染色
- ✅ 内部コントロール（他のタンパク質）正常染色

**タイムライン**: 3-5日
**コスト**: ¥80,000-150,000（複数組織、複数抗体）

---- 

### 2.4 フローサイトメトリー（細胞表面発現）

#### 目的
生細胞の細胞表面におけるPERV Envタンパク質（gp70、p15E）の発現を定量する。

#### プロトコル

**細胞準備**:
```
1. PBMC分離:
   - ヘパリン採血（10-20 ml）
   - Ficoll-Paque density gradient centrifugation
   - 300 g、30分、室温
   - バフィーコート（PBMC層）回収
   - PBSで2回洗浄

2. 細胞数調整:
   - Trypan blue染色で生存率確認（>95%）
   - 10^6細胞/mlに調整
```

**表面染色**:
```
1. Fc receptor blocking:
   - Anti-CD16/CD32抗体（1:100）
   - 4°C、10分

2. 一次抗体反応:
   - Anti-PERV gp70-FITC (1:100)
   - Anti-PERV p15E-PE (1:100)
   - 4°C、30分、暗所

3. 洗浄:
   - FACS buffer (PBS + 2% FBS + 0.1% NaN3) × 2回

4. 生死判定染色:
   - 7-AADまたはPropidium Iodide添加
   - 氷上、5分
```

**フローサイトメトリー測定**:
```
機器設定:
- Flow cytometer: BD FACSCanto II または同等品
- Events: 最低10,000細胞（生細胞ゲート内）
- Compensation: Single-color controlsで調整

ゲーティング戦略:
1. FSC-A vs SSC-A: リンパ球ゲート
2. FSC-A vs FSC-H: シングレット選択
3. 7-AAD陰性: 生細胞選択
4. FITC vs PE: gp70+/p15E+細胞定量

解析:
- FlowJo software使用
- 陽性閾値: Isotype controlの99%タイル値
- MFI (Mean Fluorescence Intensity)計算
```

**合格基準**:
- ✅ gp70陽性細胞: \<2%（バックグラウンドレベル）
- ✅ p15E陽性細胞: \<2%
- ✅ 陽性コントロール（PK-15細胞）: \>50%陽性

**タイムライン**: 4-6時間
**コスト**: ¥50,000-80,000/サンプル

---- 

## レベル3：機能レベルの検証

### 3.1 共培養アッセイ（ゴールドスタンダード）

#### 目的
ブタ細胞から産生される感染性PERVウイルスの有無を機能的に確認する。

#### プロトコル概要

**フェーズ1：ウイルス産生（ブタ細胞培養）**

```
1. ブタ細胞分離:
   - 腎臓組織またはPBMC使用
   - Collagenase消化（腎臓の場合）
   - Ficoll分離（PBMCの場合）

2. 細胞培養:
   - 培地: DMEM + 10% FBS + Pen/Strep + L-glutamine
   - マイトジェン刺激: PHA（5 μg/ml）+ IL-2（100 U/ml）
   - 培養期間: 2-4週間
   - CO2インキュベーター: 37°C、5% CO2

3. 上清回収:
   - 週2回、培地交換時に上清回収
   - 4°C、2,000 gで10分遠心（細胞除去）
   - 0.45 μmフィルター濾過
   - -80°Cで保存またはすぐに使用
```

**フェーズ2：ヒト細胞感染試験**

```
1. ヒト細胞準備:
   標的細胞（複数使用推奨）:
   - HEK-293（最も感受性高い）
   - HT1080（線維肉腫由来）
   - 初代ヒト線維芽細胞
   - ヒトPBMC

   細胞密度: 2-5 × 10^5細胞/well（24-wellプレート）

2. 感染:
   - ブタ細胞上清500 μl添加
   - Polybrene（8 μg/ml）添加（感染効率向上）
   - 24時間インキュベート

3. 洗浄:
   - PBS × 3回（未結合ウイルス除去）
   - 新鮮培地に交換

4. 長期培養:
   - 8週間継続培養
   - 週2回、培地50%交換
   - 継代: confluent時に1:3分割
```

**フェーズ3：感染検出（多重法）**

**3.1 RT-PCR（週次モニタリング）**

```python
# 週次検査プロトコル
weekly_testing = {
    "RNA抽出": {
        "sample": "細胞上清1 ml",
        "method": "QIAamp Viral RNA Mini Kit",
        "elution": "50 μl RNase-free water"
    },
    "RT-PCR": {
        "target": "PERV pol gene",
        "primers": {
            "forward": "5'-CTCCTGTTTCCTAGCCCTATT-3'",
            "reverse": "5'-GGAATTATACCTGCCATCAGC-3'"
        },
        "product_size": "285 bp",
        "cycles": 40,
        "sensitivity": "~10 copies/reaction"
    },
    "Nested PCR（確認用）": {
        "template": "1st PCR産物1 μl",
        "inner_primers": "pol遺伝子内部プライマー",
        "product_size": "187 bp",
        "sensitivity": "<5 copies/reaction"
    }
}
```

**3.2 逆転写酵素（RT）活性測定**

```
1. サンプル準備:
   - 細胞上清を100,000 gで2時間超遠心
   - ウイルスペレット回収
   - Lysis bufferで再懸濁

2. RT活性測定（ELISA法）:
   - Colorimetric RT assay kit使用
   - Poly(A)テンプレート、oligo(dT)プライマー
   - BrdU標識dUTP取り込み
   - Anti-BrdU-POD抗体で検出
   - ABTS基質で発色（405 nm吸光度測定）

3. 判定:
   - 陰性コントロール（mock培養）の2倍以上 → RT活性陽性
   - 感度: 0.01-0.1 RT units（MLV RT相当）
```

**3.3 プロウイルスDNA検出（8週後）**

```
1. ゲノムDNA抽出:
   - 感染8週後のヒト細胞から抽出
   - DNeasy Blood & Tissue Kit使用

2. PERV統合DNA検出PCR:
   - Target: PERV gag-pol junction
   - Nested PCR（高感度）
   - 1st round: 35 cycles
   - 2nd round: 30 cycles

3. ブタDNA混入除外:
   - Sus scrofa cytochrome b遺伝子PCR
   - 陰性確認（ブタ細胞混入なし）
```

**陽性/陰性コントロール**

```
必須コントロール:
1. 陽性コントロール:
   - PK-15細胞上清（既知PERV産生）
   - PERV疑似粒子（VLP）
   - → 全てのアッセイで検出されること

2. 陰性コントロール:
   - Mock感染（培地のみ）
   - 未感染ヒト細胞
   - → 全てのアッセイで陰性

3. システムコントロール:
   - ヒトHousekeeping gene（GAPDH）PCR陽性確認
   - 細胞生存率 >90%確認
```

**合格基準（PERV-KOブタ）**

```
8週間の観察期間において:
- ✅ RT-PCR: 全週で陰性（Ct値なしまたは>40）
- ✅ RT活性: バックグラウンドレベル（<0.01 units）
- ✅ プロウイルスDNA: PCR陰性
- ✅ ブタDNA混入: 陰性（交差汚染なし）
- ✅ 陽性コントロール: 全アッセイで明確な陽性
```

**タイムライン**: 10-12週間（ブタ細胞培養2-4週 + ヒト細胞観察8週）
**コスト**: ¥200,000/動物

---- 

### 3.2 ヒト細胞株感染パネル

#### 目的
複数のヒト細胞種を用いてPERV感染性を包括的に評価する。

#### プロトコル

**細胞株パネル**:
```
1. HEK-293（ヒト胎児腎臓）:
   - PERV-Aに最も感受性高い
   - Reference cell line

2. HT1080（ヒト線維肉腫）:
   - PERV研究の標準細胞株
   - PERV-A、PERV-B両方に感受性

3. TE671（ヒト横紋筋肉腫）:
   - 中等度感受性

4. 初代ヒト線維芽細胞:
   - 生理的関連性高い
   - ドナー由来変動あり

5. ヒトPBMC:
   - 臨床的に最も重要
   - PHA/IL-2刺激必要
```

**並行感染試験**:
```
同一ブタ細胞上清で5種類のヒト細胞を同時感染
→ 各細胞株で8週間観察
→ RT-PCR、RT活性、プロウイルスDNA検出

解析:
- 細胞株間の感受性差を評価
- 最も感受性の高い細胞株でも陰性 → PERV-free確定
```

**タイムライン**: 10-12週間
**コスト**: ¥400,000-600,000（5細胞株並行）

---- 

## 細胞表面タンパク質の発現メカニズム

### PERVライフサイクルと表面タンパク質

```
【レトロウイルス複製サイクル】

1. ウイルス侵入
   ↓
2. 逆転写（RNA → DNA）
   ↓
3. 核移行・ゲノム統合
   ↓
4. 転写（プロウイルスDNA → mRNA）
   ↓
5. 翻訳
   ├─ Gag前駆体（Pr65）→ p27, p10（細胞質）
   ├─ Pol前駆体（PR, RT, IN）（細胞質）
   └─ Env前駆体（gPr80）→ gp70, p15E（小胞体）
   ↓
6. Envタンパク質プロセシング
   ├─ 小胞体でシグナルペプチド切断
   ├─ gp70（SU）とp15E（TM）にプロセシング
   ├─ N-グリコシル化
   └─ ゴルジ体を経て細胞膜へ輸送
   ↓
7. 細胞膜にgp70/p15E発現 ← ここが検出ポイント
   ↓
8. ウイルス粒子出芽（Budding）
   ├─ Gagポリタンパク質が細胞膜下に集合
   ├─ 細胞膜のgp70/p15Eを取り込み
   └─ 出芽によりエンベロープウイルス形成
   ↓
9. 成熟（Maturation）
   └─ プロテアーゼによりGagが切断 → p27, p10
```

### なぜ細胞表面にタンパク質が発現するのか

**エンベロープウイルスの特性**:
```
PERVはエンベロープウイルス（膜を持つ）
→ 宿主細胞の細胞膜を利用して出芽
→ 出芽に必要なEnvタンパク質（gp70, p15E）が細胞表面に発現

非エンベロープウイルス（例: AAV）:
→ 細胞溶解で放出
→ 細胞表面タンパク質発現なし
```

**感染細胞の特徴**:
```
PERV感染細胞:
├─ 細胞表面にgp70（受容体結合ドメイン）
├─ 細胞表面にp15E（膜貫通ドメイン）
├─ これらがウイルス粒子に組み込まれる
└─ フローサイトメトリー、免疫蛍光で検出可能

PERV-KO細胞:
├─ Env遺伝子が破壊
├─ gp70/p15Eタンパク質が合成されない
├─ 細胞表面に発現なし
└─ ウイルス粒子形成不能（たとえGag/Polがあっても）
```

### 検出方法の原理

**フローサイトメトリー**:
```python
detection_principle = {
    "標的": "生細胞の細胞表面gp70/p15E",
    "抗体": "蛍光標識抗PERV抗体",
    "結合": "細胞表面のEnvタンパク質に特異的結合",
    "検出": "蛍光シグナル強度を個別細胞レベルで定量",
    "結果": {
        "野生型ブタPBMC": "30-70%の細胞がgp70陽性",
        "PERV-KOブタPBMC": "<2%（バックグラウンド）",
        "感染ヒト細胞": "10-50%陽性（感染成立時）"
    }
}
```

**免疫蛍光顕微鏡**:
```
1. 生細胞を抗PERV抗体で染色（非透過処理）
2. 細胞表面のみ染色される
3. 共焦点顕微鏡で細胞膜の蛍光観察
4. 野生型: 細胞膜全体に蛍光リング
5. PERV-KO: 蛍光なし
```

---- 

## 多世代安定性検証

### 目的
CRISPR編集が生殖細胞系列を通じて安定に遺伝することを確認する。

### プロトコル

**世代間比較**:
```
F0世代（創始世代）:
└─ Tier 3包括的検証（¥730,000）
   ├─ WGS（全ゲノム）
   ├─ ddPCR
   ├─ Western Blot（複数組織）
   ├─ 共培養アッセイ
   └─ フローサイトメトリー

F1世代（第一子孫）:
└─ Tier 2標準検証（¥430,000）
   ├─ WGS（遺伝継承確認）
   ├─ ddPCR（コピー数）
   ├─ Western Blot
   └─ 共培養アッセイ

F2世代（第二子孫、臨床使用）:
└─ Tier 2標準検証（¥430,000）
   ├─ WGS
   ├─ ddPCR
   ├─ Western Blot
   └─ 共培養アッセイ
```

**遺伝パターン解析**:
```python
# メンデル遺伝の確認
inheritance_analysis = {
    "F0 × F0 (両方ホモKO)": {
        "expected": "F1全てホモKO",
        "validation": "WGSで全F1個体のPERV遺伝子座確認"
    },
    "F1 × F1 (両方ホモKO)": {
        "expected": "F2全てホモKO",
        "validation": "ddPCRで0コピー確認"
    },
    "世代間一貫性": {
        "DNA level": "全世代で同一の編集パターン",
        "Protein level": "全世代でPERVタンパク質不在",
        "Functional level": "全世代で共培養陰性"
    }
}
```

**合格基準**:
- ✅ F0、F1、F2で一貫したPERV-KO表現型
- ✅ 予期しない表現型変異なし
- ✅ 予期しない遺伝パターンなし（モザイク、復帰変異など）

---- 

## 品質管理とドキュメンテーション

### 必須記録項目

```
各検証実施時の記録:
1. 実施日
2. 実施者
3. サンプル情報（個体ID、組織種類、採取日）
4. 試薬ロット番号（全試薬）
5. 機器情報（較正日、メンテナンス記録）
6. 陽性/陰性コントロール結果
7. 生データ（画像、クロマトグラム、配列データ）
8. 解析パラメータ
9. 逸脱記録（プロトコルからの逸脱があれば）
10. 結果判定と判定者署名
```

### PMDA提出用データパッケージ

```
提出書類構成:
├─ 検証サマリーレポート
├─ 個別検証結果
│  ├─ WGS報告書（VCFファイル、アライメント統計）
│  ├─ ddPCR定量データ（生データ、校正曲線）
│  ├─ Western Blot画像（生画像、定量データ）
│  ├─ 共培養アッセイ結果（週次PCRデータ、8週サマリー）
│  └─ フローサイトメトリーデータ（FCSファイル、ゲーティング戦略）
├─ プロトコル文書（SOP）
├─ 試薬・機器リスト
├─ 品質管理記録
│  ├─ 陽性/陰性コントロール結果
│  ├─ 機器較正記録
│  └─ 試薬バリデーション
└─ 逸脱報告（あれば）
```

---- 

## トラブルシューティング

### Western Blotでバンドが見えない（PERV-KOブタ以外で）

**原因と対策**:
```
1. タンパク質濃度不足:
   → BCA assayで再定量、50 μg/レーンまで増量

2. 抗体の問題:
   → 陽性コントロール（PK-15）で抗体機能確認
   → 別ロットまたは別ベンダーの抗体試行

3. 転写効率低下:
   → Ponceau S染色で転写確認
   → Transfer条件最適化（電圧、時間）

4. ブロッキング過剰:
   → BSA濃度を5% → 3%に低下
   → ブロッキング時間短縮

5. ECL感度不足:
   → 高感度ECL substrate使用
   → 露光時間延長（最大30分）
```

### 共培養アッセイで陽性コントロールが陰性

```
1. PK-15細胞の問題:
   → 継代数確認（高継代でPERV産生低下）
   → 新しいストック解凍

2. ヒト細胞の感受性低下:
   → HEK-293の継代数確認
   → 新しいロット入手

3. Polybrene濃度:
   → 8 μg/mlに調整（高すぎると細胞毒性）

4. ウイルス力価不足:
   → PK-15培養期間延長（2 → 4週間）
   → マイトジェン刺激追加
```

### ddPCRでドロップレット数不足

```
1. ドロップレット生成失敗:
   → Oil/sampleの量確認
   → ガスケット交換
   → Droplet Generatorクリーニング

2. サンプルDNA濃度不適切:
   → 推奨範囲: 10-100 ng/反応
   → 希釈系列で最適化

3. プライマー/プローブの問題:
   → 標準PCRで増幅確認
   → 再設計またはIDT注文確認
```

---- 

## まとめ

### 3レベル検証の重要性

```
レベル1（DNA）: 遺伝子が存在しない
    +
レベル2（タンパク質）: タンパク質が発現していない
    +
レベル3（機能）: 感染性ウイルスが産生されない
    ↓
完全なPERVフリーの証明
```

### 最小限vs推奨vs包括的

| 検証レベル     | 最小限（Tier 1） | 推奨（Tier 2）   | 包括的（Tier 3）          |
| --------- | ----------- | ------------ | -------------------- |
| **DNA**   | WGS + ddPCR | WGS + ddPCR  | WGS + ddPCR + Sanger |
| **タンパク質** | なし          | Western Blot | Western + IHC + Flow |
| **機能**    | 共培養         | 共培養          | 共培養 + 細胞株パネル         |
| **コスト**   | ¥380,000    | ¥430,000     | ¥730,000             |
| **規制リスク** | 10-15%      | \\\<5%       | 0%                   |

### 次のステップ

これらの検証方法を実装するには:
1. [コストパフォーマンス分析][1]で予算計画
2. [規制要件][2]でコンプライアンス確認
3. [ELISA分析][3]で代替方法検討

---- 

## クローニング個体の特殊な検証プロトコル

### クローニングと遺伝的ドリフトのモニタリング

#### 目的
体細胞クローニングで生産されたPERV-KOブタの遺伝的同一性と安定性を確認する。

### クローニング技術バリデーション（初期20頭）

#### プロトコル

**Phase 1: Founder（eGenesis）との比較解析**

```python
# 全20頭にWGS実施後の解析
cloning_validation = {
    "variant_calling": {
        "method": "GATK HaplotypeCaller",
        "reference": "Founder pig genome (from eGenesis)",
        "comparison": "Each clone vs Founder"
    },

    "metrics": {
        "SNP_concordance": ">99.99% 期待",
        "indel_concordance": ">99.99% 期待",
        "PERV_loci_identity": "100% 必須",
        "off_target_sites": "Founderと完全一致"
    },

    "acceptable_differences": {
        "total_SNPs": "<10 個/ゲノム",
        "PERV_region_SNPs": "0個（許容なし）",
        "functional_gene_variants": "0個（許容なし）",
        "mitochondrial_variants": "許容（核ゲノム由来でない）"
    }
}
```

**Phase 2: クローン間の一貫性評価**

```python
# 20頭のクローン間比較
inter_clone_analysis = {
    "pairwise_comparison": {
        "method": "All-vs-all SNP calling",
        "expected": "99.99%+ identity between clones",
        "threshold": "Discordant sites <20 per pair"
    },

    "clustering_analysis": {
        "method": "PCA on genetic variants",
        "expected": "Tight clustering of all 20 clones",
        "outlier_detection": "Any clone >2 SD from mean"
    },

    "PERV_specific_check": {
        "method": "Multiple sequence alignment of 62 PERV loci",
        "expected": "100% identity across all clones",
        "action_if_discordant": "Exclude clone, investigate cloning process"
    }
}
```

**合格基準（技術バリデーション）**:
- ✅ 全20頭がFounderと\>99.99%一致
- ✅ PERV遺伝子座は100%一致（全クローン）
- ✅ Off-target部位に予期しない変異なし
- ✅ クローン間の遺伝的ばらつき \<0.01%

**不合格時の対応**:
- 不一致率\>0.01%のクローン: 除外、追加クローン作製
- PERV遺伝子座の不一致: プロセス全体の再検証
- Off-target変異の出現: クローニング手技の見直し

---- 

### 代表サンプリング検証（生産拡大フェーズ）

#### ロット管理とサンプリング戦略

```python
lot_based_sampling = {
    "lot_definition": {
        "criteria": "同一クローニングバッチ（1-2週間以内）",
        "typical_size": "10-20頭",
        "documentation": "ロット番号、作製日、使用細胞ライン記録"
    },

    "sampling_frequency": {
        "early_production": {
            "lots": "1-10（頭数21-100）",
            "sampling": "各ロットから1頭（10%）",
            "rationale": "クローニング安定性の継続確認"
        },
        "established_production": {
            "lots": "11以降（頭数101+）",
            "sampling": "20頭ごとに1頭（5%）",
            "condition": "初期データで忠実度実証済み",
            "exception": "異常検出時は全頭検査に復帰"
        }
    },

    "wgs_analysis": {
        "comparison": "Founder + 初期20頭の平均",
        "metrics": "SNP concordance, PERV loci identity",
        "trending": "累積SNP数の経時変化プロット",
        "alert_threshold": "Cumulative SNPs >50"
    }
}
```

#### 全頭簡易検証プロトコル

**ddPCR（全クローン頭実施）**:

```python
# クローニング特化のddPCRプロトコル
clone_ddpcr = {
    "purpose": "遺伝的同一性のスクリーニング",

    "target_genes": {
        "PERV_pol": {
            "expected": "0 copies/genome",
            "action_if_positive": "即WGS実施、クローニング異常調査"
        },
        "beta_actin": {
            "expected": "2 copies/genome (diploid)",
            "QC": "1.8-2.2 範囲外は再測定"
        }
    },

    "additional_markers": {
        "PERV_A_env": "0 copies",
        "PERV_B_env": "0 copies",
        "PERV_C_env": "0 copies",
        "rationale": "サブタイプ特異的確認"
    },

    "clone_specific_threshold": {
        "LOD": "<0.1 copies/genome",
        "stricter_than_standard": "通常0.5の5倍感度",
        "reason": "クローニング変異の早期検出"
    }
}
```

**マルチプレックスPCR（PERV-A/B/C同時検出）**:

```
プライマー設計:
- PERV-A env (RBD): 200 bp product
- PERV-B env (RBD): 250 bp product
- PERV-C env (RBD): 300 bp product
- Sus scrofa GAPDH: 150 bp (内部コントロール)

反応条件:
- アニーリング温度: 58°C
- 増幅サイクル: 35 cycles
- 検出: アガロースゲル電気泳動

結果解釈:
- 全クローン: GAPDH陽性、PERV-A/B/C全て陰性
- PERV陽性の場合: 直ちにWGS + クローニング記録確認
```

---- 

### クローニング異常の検出と対応

#### 異常パターンの同定

```python
anomaly_detection = {
    "genetic_anomalies": {
        "unexpected_SNPs": {
            "definition": "Founderと>10 SNP差",
            "frequency_threshold": "2%以上の頻度で発生",
            "action": [
                "WGS full analysis",
                "クローニング手技レビュー",
                "使用細胞ラインの品質確認"
            ]
        },

        "PERV_reversion": {
            "definition": "PERVコピー数>0",
            "probability": "<0.001% (理論的)",
            "action": [
                "即座に該当クローン除外",
                "全ロット緊急WGS",
                "eGenesisへの報告",
                "PMDA緊急報告（24時間以内）"
            ]
        },

        "chromosomal_abnormality": {
            "detection": "WGSカバレッジ異常、大規模欠失/重複",
            "action": "除外、クローニングプロトコル見直し"
        }
    },

    "phenotypic_anomalies": {
        "growth_retardation": "体重<正常範囲の80%",
        "organ_malformation": "画像診断で検出",
        "action": "WGS + 臨床検査、使用中止"
    },

    "epigenetic_concerns": {
        "note": "WGSでは検出不可",
        "monitoring": "遺伝子発現プロファイル（任意）",
        "functional_assay": "共培養アッセイが最終確認"
    }
}
```

#### 段階的エスカレーション

```
レベル1（低リスク異常）:
- 1-2頭でddPCR異常値（β-actin範囲外など）
- 対応: 再測定、個体追跡、次回ロットで注意

レベル2（中リスク異常）:
- 1ロット内で2頭以上の異常
- 対応: ロット全体の代表サンプルWGS（3頭）、次ロット全頭検査

レベル3（高リスク異常）:
- PERVコピー数陽性
- 複数ロットで遺伝的不一致
- 対応:
  1. 全生産停止
  2. 既存全クローンの緊急WGS
  3. クローニング技術の完全再バリデーション
  4. PMDA緊急報告
```

---- 

### クローニング検証のPMDA報告書フォーマット

#### 必須セクション

```markdown
# 体細胞クローニング技術バリデーション報告書

## 1. エグゼクティブサマリー
- Founder pig由来、初期20頭クローンのWGS完了
- 遺伝的忠実度: 99.99%以上（平均SNP差: 3.2個/ゲノム）
- PERV遺伝子座: 100%一致（全20頭×62遺伝子座）
- Off-target部位: Founderと完全一致
- 結論: クローニング技術は遺伝的に安定

## 2. 技術詳細
### 2.1 クローニング方法
- 体細胞核移植（SCNT）
- ドナー細胞: Founder pig線維芽細胞
- レシピエント卵子: 野生型ブタ
- 施設: [クローニング実施施設名]
- 成功率: XX% (Y頭/Z試行)

### 2.2 WGS解析パラメータ
- カバレッジ: 30× (範囲: 28-35×)
- マッピング率: >98%
- バリアントコール: GATK 4.2
- 参照ゲノム: Sscrofa11.1 + Founder WGS

## 3. 遺伝的忠実度データ
### 3.1 SNP解析
[表: 各クローンのSNP数、PERV領域のSNP数、Off-target変異数]

### 3.2 クローン間一貫性
[PCAプロット: 20クローンの遺伝的クラスタリング]
[ヒートマップ: Pairwise genetic distance]

## 4. 代表サンプリング戦略の正当性
### 4.1 統計学的根拠
- サンプルサイズ計算: 10%サンプリングで95%信頼区間
- 検出力: 0.1%の変異を95%確率で検出

### 4.2 継続的モニタリング計画
- 10頭ごとに1頭WGS（生産フェーズ）
- 全頭ddPCR（PERVコピー数）
- 異常検出時の全頭検査プロトコル

## 5. 品質保証
### 5.1 逸脱管理
[クローニング失敗、異常個体の記録と対応]

### 5.2 トレーサビリティ
[個体ID → クローニング日 → 使用細胞ライン → Founder]

## 6. 結論
体細胞クローニング技術は99.99%以上の遺伝的忠実度を実証。
代表サンプリング戦略は科学的に妥当であり、全頭WGSは不要と判断。
```

---- 

### クローニング vs 繁殖の比較

| 項目            | 体細胞クローニング       | 自然繁殖            | 推奨         |
| ------------- | --------------- | --------------- | ---------- |
| **遺伝的均一性**    | 極めて高い（\>99.99%） | 低い（メンデル分離）      | クローニング     |
| **PERV遺伝子座**  | 100%一致          | 50%遺伝（ヘテロ接合の場合） | クローニング     |
| **WGS頻度**     | 代表サンプリング可       | 全頭必要            | クローニング     |
| **検証コスト**     | 低（全頭簡易検証）       | 高（全頭完全検証）       | クローニング     |
| **生産スピード**    | 速い              | 遅い（妊娠・成熟期間）     | クローニング     |
| **エピジェネティック** | リスクあり（要監視）      | 正常              | 繁殖         |
| **総合評価**      | 異種移植用途に最適       | 研究用途            | **クローニング** |


[1]:	PERV_COST_ANALYSIS_JP.md
[2]:	PERV_REGULATORY_REQUIREMENTS_JP.md
[3]:	PERV_ELISA_ANALYSIS_JP.md