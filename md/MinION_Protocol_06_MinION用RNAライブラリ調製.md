# MinION用メタゲノム解析プロトコル
# 第6章: MinION用RNAライブラリ調製

## 目次

1. [RNAライブラリ調製概要](#1-rnaライブラリ調製概要)
2. [RNA病原体とcDNA変換の重要性](#2-rna病原体とcdna変換の重要性)
3. [Direct RNA Sequencing vs cDNA Sequencing](#3-direct-rna-sequencing-vs-cdna-sequencing)
4. [cDNA合成プロトコル (ストランド特異的)](#4-cdna合成プロトコルストランド特異的)
5. [Second strand合成とdsDNA調製](#5-second-strand合成とdsdna調製)
6. [MinION RNAライブラリ調製 (SQK-RNA002)](#6-minion-rnaライブラリ調製-sqk-rna002)
7. [cDNA-Seq用ライブラリ調製 (SQK-LSK114)](#7-cdna-seq用ライブラリ調製-sqk-lsk114)
8. [RNA品質評価とQC基準](#8-rna品質評価とqc基準)
9. [RNA特異的トラブルシューティング](#9-rna特異的トラブルシューティング)
10. [ALCOA+準拠記録](#10-alcoa準拠記録)

---

## 1. RNAライブラリ調製概要

### 1.1 RNA病原体検出の重要性

**PMDA 91病原体リストにおけるRNA病原体**

```yaml
RNA病原体の内訳:
  RNA単独病原体: 20種
    - RNAウイルス (18種):
      * インフルエンザウイルス
      * 日本脳炎ウイルス
      * 狂犬病ウイルス
      * 豚コレラウイルス
      * 豚繁殖・呼吸障害症候群ウイルス (PRRSV)
      * エンテロウイルス
      * その他12種

  DNA + RNA両方検出可能: 11種
    - レトロウイルス (PERV含む)
    - 一部のDNAウイルス (転写RNA検出)

  合計RNA検出対象: 31種 (91種中34%)

重要性:
  - RNA病原体は活発な感染を示す指標
  - DNA病原体でも転写RNA検出で活性評価可能
  - cfRNAは感染初期段階で検出感度高い
```

### 1.2 cfRNA解析の技術的課題

**RNAの脆弱性**

```yaml
RNA固有の問題点:

  1. 不安定性:
     - RNase ubiquity (どこにでも存在)
     - 半減期: 数分〜数時間
     - 室温で急速分解

  2. 低収量:
     - 血漿中cfRNA濃度: DNA の 1/10〜1/100
     - 期待収量: 5-50 ng (cfDNAの10-100 ng比)
     - ロス: 各精製ステップで20-30%

  3. 品質管理困難:
     - RIN (RNA Integrity Number) 維持困難
     - 部分分解RNAの混在
     - 定量の不正確性

対策の必要性:
  - RNase-free環境の徹底
  - 迅速な処理 (抽出からライブラリ調製まで)
  - 高感度ライブラリ調製キット使用
  - 適切な品質評価指標
```

### 1.3 MinION RNA-Seqの2つのアプローチ

**戦略比較**

```yaml
アプローチ1: Direct RNA Sequencing
  使用キット: SQK-RNA002

  利点:
    - RNAを直接シーケンス (cDNA変換不要)
    - ネイティブRNA修飾検出可能
    - ストランドバイアスなし

  欠点:
    - 必要RNA量: 500 ng (高い)
    - スループット: DNAの1/5程度
    - Poly(A) tail必須 (mRNA向け)
    - ウイルスRNA検出に不向き

  PMDA病原体検出評価: △ (RNA量要求高すぎ)

アプローチ2: cDNA Sequencing
  使用キット:
    - cDNA合成: NEBNext Ultra II Directional RNA Library Prep
    - MinIONライブラリ: SQK-LSK114

  利点:
    - 低RNA量対応: 10-100 ng
    - 高スループット (DNAと同等)
    - 全RNAタイプ対応 (mRNA, rRNA, viral RNA)
    - 安定 (cDNA化後の分解なし)

  欠点:
    - RNA修飾情報喪失
    - cDNA合成ステップ追加 (時間・コスト)
    - PCR biasの可能性

  PMDA病原体検出評価: ◎ (推奨)

推奨戦略:
  cfRNA → cDNA変換 → MinION LSK114ライブラリ調製
```

---

## 2. RNA病原体とcDNA変換の重要性

### 2.1 対象RNA病原体の分類

**RNAウイルス分類 (ゲノムタイプ別)**

```yaml
(+)鎖RNAウイルス (Positive-sense RNA):
  特徴: ゲノムRNAがそのままmRNAとして機能

  対象病原体:
    - 豚コレラウイルス (Classical swine fever virus)
      ゲノム: 12.3 kb, (+)ssRNA
      科: Flaviviridae

    - 豚繁殖・呼吸障害症候群ウイルス (PRRSV)
      ゲノム: 15.4 kb, (+)ssRNA
      科: Arteriviridae

    - エンテロウイルス (Enterovirus)
      ゲノム: 7.5 kb, (+)ssRNA
      科: Picornaviridae

    - 日本脳炎ウイルス (Japanese encephalitis virus)
      ゲノム: 11 kb, (+)ssRNA
      科: Flaviviridae

    - 豚流行性下痢ウイルス (PEDV)
      ゲノム: 28 kb, (+)ssRNA
      科: Coronaviridae

  cDNA合成での注意:
    - ゲノムRNAを直接テンプレートに可能
    - Random hexamer primerで効率的
    - Poly(A) tailなしでも合成可能

(-)鎖RNAウイルス (Negative-sense RNA):
  特徴: 相補鎖がmRNA、ゲノムは鋳型鎖

  対象病原体:
    - インフルエンザAウイルス (Influenza A virus)
      ゲノム: 8セグメント, 合計13.5 kb, (-)ssRNA
      科: Orthomyxoviridae

    - 狂犬病ウイルス (Rabies virus)
      ゲノム: 12 kb, (-)ssRNA
      科: Rhabdoviridae

    - センダイウイルス (Sendai virus)
      ゲノム: 15.4 kb, (-)ssRNA
      科: Paramyxoviridae

  cDNA合成での注意:
    - 相補鎖とゲノム鎖の両方検出必要
    - ストランド特異的プライマー推奨
    - Viral replicationで(+)鎖mRNA産生

2本鎖RNAウイルス (dsRNA):
  対象病原体:
    - ロタウイルス (Rotavirus)
      ゲノム: 11セグメント, 合計18.6 kb, dsRNA
      科: Reoviridae

  cDNA合成での注意:
    - dsRNA変性処理必要
    - Random primerで両鎖増幅

レトロウイルス (RNA → DNA逆転写):
  対象病原体:
    - 豚内因性レトロウイルス (PERV)
      ゲノム: 9 kb, (+)ssRNA (逆転写酵素保有)
      科: Retroviridae

  cDNA合成での注意:
    - ゲノムDNA統合型も存在
    - RNA検出で活性評価
    - cfDNAとcfRNA両方で検出
```

### 2.2 ストランド特異性の重要性

**ストランド特異的シーケンスの必要性**

```yaml
非ストランド特異的cDNA合成の問題:
  - (+)鎖と(-)鎖の区別不可
  - ゲノム鎖 vs 転写鎖の識別不能
  - 遺伝子発現方向が不明
  - ウイルス複製ステージ判定困難

ストランド特異的cDNA合成の利点:
  1. 正確なRNA定量:
     - (+)鎖RNA量と(-)鎖RNA量を別々に定量
     - ウイルス複製活性の評価

  2. 転写方向の同定:
     - 遺伝子アノテーション精度向上
     - アンチセンス転写検出

  3. RNA編集・修飾の検出:
     - 方向特異的変異解析
     - エピトランスクリプトーム解析

PMDA病原体検出における重要性:
  - 活発な感染 vs 潜伏感染の区別
  - ワクチン株 vs 野生株の識別 (転写パターン差)
  - 薬剤耐性株の検出 (特定遺伝子発現)
```

---

## 3. Direct RNA Sequencing vs cDNA Sequencing

### 3.1 詳細比較表

| 項目 | Direct RNA-Seq | cDNA-Seq |
|-----|---------------|----------|
| **必要RNA量** | 500 ng | 10-100 ng |
| **cfRNA対応** | ✗ (量不足) | ◎ |
| **Poly(A) tail要求** | 必須 | 不要 |
| **ウイルスRNA検出** | △ (mRNAウイルスのみ) | ◎ (全RNA) |
| **RNA修飾検出** | ◎ (m6A, pseudouridineなど) | ✗ |
| **ストランド特異性** | ◎ (ネイティブ) | ○ (プロトコル依存) |
| **スループット** | 0.5-2 Gb/flowcell | 15-30 Gb/flowcell |
| **リード長** | 500-3,000 nt (RNA依存) | 5-50 kb (cDNA依存) |
| **精度** | Q10-15 | Q20-30 (Duplex) |
| **ライブラリ調製時間** | 3時間 | 8-12時間 |
| **コスト** | ¥50,000/sample | ¥40,000/sample |
| **PMDA適用性** | ✗ | ◎ |

### 3.2 cfRNA血漿サンプルでの推奨戦略

**結論: cDNA-Seqアプローチ一択**

```yaml
理由:
  1. RNA量制限:
     cfRNA収量: 5-50 ng (平均20 ng)
     Direct RNA要求: 500 ng
     → 10-100倍不足、実施不可能

  2. Poly(A) tail非保有RNA:
     ウイルスRNA多くはPoly(A) tailなし
     → Direct RNAではキャプチャ不可

  3. スループット要求:
     91病原体検出には深いカバレッジ必要
     → cDNA-Seqの高スループット必須

  4. 安定性:
     cDNA化後は分解リスクなし
     → 長期保管、再解析可能

推奨プロトコル:
  cfRNA抽出 (Zymo Kit)
    ↓
  DNase処理 (DNA残存除去)
    ↓
  ストランド特異的cDNA合成 (NEBNext)
    ↓
  Second strand合成
    ↓
  MinION LSK114ライブラリ調製
    ↓
  シーケンス (48時間、Duplex mode)
```

---

## 4. cDNA合成プロトコル (ストランド特異的)

### 4.1 必要試薬・機器

**試薬キット**

```yaml
推奨キット: NEBNext Ultra II Directional RNA Library Prep Kit
  カタログ番号: E7760S/L
  反応数: 24 reactions (S) / 96 reactions (L)

  キット内容:
    - NEBNext RNA First Strand Synthesis Buffer (10X)
    - NEBNext Random Primers
    - NEBNext Strand Specificity Reagent (SSR)
    - ProtoScript II Reverse Transcriptase
    - NEBNext Second Strand Synthesis Reaction Buffer (10X)
    - NEBNext Second Strand Synthesis Enzyme Mix

  価格: ¥85,000 (24 reactions)

追加試薬:
  - TURBO DNase (Thermo Fisher, AM2238)
    → RNA中の残存DNA完全除去
    価格: ¥45,000 (50 reactions)

  - RNaseOUT Recombinant RNase Inhibitor (Invitrogen, 10777019)
    → RNA分解防止
    価格: ¥35,000 (5,000 units)

  - Nuclease-free water
  - 100% エタノール (分子生物学グレード)
  - AMPure XP Beads (Beckman Coulter, A63881)

機器:
  - サーマルサイクラー (正確な温度制御)
  - マグネティックスタンド (1.5 mL tube用)
  - ボルテックスミキサー
  - マイクロ遠心機
  - 氷・ドライアイス
```

### 4.2 DNase処理プロトコル

**目的: cfRNA中の残存DNA完全除去**

RNAとDNAは同時抽出されるため、RNA特異的解析にはDNA除去が必須です。

```yaml
TURBO DNase処理プロトコル:

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ1: DNase反応系の調製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  反応系 (20 μL):
    cfRNA (Zymo Kit溶出液): 15 μL
    TURBO DNase Buffer (10X): 2 μL
    TURBO DNase (2 U/μL): 2 μL
    RNaseOUT (40 U/μL): 1 μL
    ────────────────────────
    総量: 20 μL

  手順:
    1. 氷上で全試薬を調製
    2. 0.2 mL PCRチューブに以下の順で添加:
       a. cfRNA: 15 μL
       b. TURBO DNase Buffer (10X): 2 μL
       c. RNaseOUT: 1 μL
       d. TURBO DNase: 2 μL (最後に添加)

    3. ピペッティング混合: 10回 (優しく)
    4. スピンダウン: 3,000 rpm、3秒

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ2: DNase反応
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  インキュベーション:
    温度: 37°C
    時間: 30分
    機器: サーマルサイクラー

  重要:
    - 蓋加熱: 50°C (蒸発防止)
    - 反応終了まで静置

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ3: DNase不活化と精製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  不活化試薬添加:
    DNase Inactivation Reagent: 2 μL (kit付属)
    混合: ボルテックス 5秒
    インキュベート: 室温、5分 (時々タッピング)
    遠心: 10,000×g、1.5分

  上清回収:
    - 新しい1.5 mL LoBind tubeへ上清18 μLを移す
    - 沈殿は触れない (DNaseとDNA複合体)

  または AMPure XP精製 (推奨):
    1. AMPure XP beads添加: 36 μL (2.0× ratio)
    2. ピペッティング混合: 10回
    3. 室温インキュベート: 5分
    4. マグネティックスタンド: 5分
    5. 上清除去
    6. 80% エタノール洗浄: 2回 (各200 μL)
    7. 風乾: 5分 (過乾燥注意)
    8. Nuclease-free water溶出: 11 μL
    9. 室温インキュベート: 2分
    10. マグネティックスタンド: 2分
    11. 上清10 μL回収 → 新チューブ

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ4: DNase処理効果確認 (QC)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Qubit測定:
    - Qubit RNA HS Assay: RNA濃度測定
    - Qubit DNA HS Assay: 残存DNA測定

    合格基準:
      RNA回収率: >80% (処理前比)
      残存DNA: <1% (処理前比)
      DNA/RNA比: <0.01

  qPCR検証 (オプション、推奨):
    豚ゲノムDNA特異的プライマー使用
    (例: β-actin gene)

    反応系:
      Template: DNase処理RNA 1 μL
      Control: 処理前RNA 1 μL (1/100希釈)
      Primer: Sus scrofa β-actin

    合格基準:
      ΔCt >10 (処理前比で>1,000倍減少)

  記録:
    - 処理前RNA濃度: _____ ng/μL
    - 処理後RNA濃度: _____ ng/μL
    - 回収率: _____ %
    - 残存DNA: _____ ng/μL
    - DNA/RNA比: _____
    - 判定: □合格 □不合格
```

**重要な注意事項**

```yaml
RNase汚染防止:
  - 全ステップで新しいフィルターチップ使用
  - RNaseZAP処理済みピペット使用
  - 氷は専用のRNase-free ice作成
  - 手袋は頻繁に交換

迅速な処理:
  - DNase処理からcDNA合成まで連続実施
  - 中断する場合: -80°C保存、1週間以内

品質管理:
  - Bioanalyzer RIN測定 (理想: RIN >7)
  - 分解が進んでいる場合でもcDNA合成可能
  - ただし、収量とカバレッジが低下
```

### 4.3 First Strand cDNA合成

**ストランド特異的cDNA合成の原理**

NEBNext Ultra II Directional キットは、dUTP法を使用してストランド特異性を実現します:

```yaml
原理:
  1. First strand合成:
     RNA → cDNA (Random primerまたはOligo(dT))
     dNTP使用 (通常のdATP, dCTP, dGTP, dTTP)

  2. RNA分解:
     RNase H処理でテンプレートRNA除去

  3. Second strand合成:
     dUTP含有dNTP mix使用
     → Second strandにdUTPが取り込まれる

  4. Second strand選択的分解:
     UDG (Uracil-DNA Glycosylase)処理
     → dUTP含有Second strandのみ除去
     → First strand (元のRNA情報保持) のみ残存

結果:
  元のRNA方向を保持したcDNAライブラリ
  (+)鎖RNAと(-)鎖RNAを区別可能
```

**First Strand合成プロトコル**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ1: RNA fragmentation (オプション)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  目的: 長鎖RNAの断片化 (>2 kb RNA)

  ※ cfRNAの場合、通常すでに断片化されているため省略可
  ※ ウイルスゲノムRNA (>10 kb) がある場合は実施推奨

  条件 (実施する場合):
    温度: 94°C
    時間: 15分 (目標500-2,000 nt fragment)
    バッファー: NEBNext First Strand Synthesis Buffer

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ2: First Strand反応系調製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  反応系 (20 μL):
    DNase処理RNA: 8 μL (20-100 ng推奨)
    NEBNext Random Primers: 2 μL
    NEBNext First Strand Synthesis Buffer (10X): 2 μL
    NEBNext Strand Specificity Reagent (SSR): 8 μL
      (Actinomycin D含有、rRNA由来cDNA抑制)
    ────────────────────────────
    総量: 20 μL

  調製手順:
    1. 氷上で0.2 mL PCRチューブに添加:
       a. DNase処理RNA: 8 μL
       b. NEBNext Random Primers: 2 μL

    2. プライマーアニーリング:
       サーマルサイクラー設定:
         65°C - 5分
         氷上 - 1分以上

    3. 氷上で追加試薬添加:
       c. First Strand Synthesis Buffer (10X): 2 μL
       d. Strand Specificity Reagent (SSR): 8 μL
       ※ SSRは粘性高い、ゆっくりピペッティング

    4. ピペッティング混合: 10回
    5. スピンダウン: 3,000 rpm、5秒

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ3: cDNA合成反応
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  サーマルサイクラー設定:
    ステップ1: 25°C - 10分 (Random primer extension開始)
    ステップ2: 42°C - 15分 (cDNA合成)
    ステップ3: 70°C - 15分 (酵素失活)
    ステップ4: 4°C - ホールド

  重要:
    - 蓋加熱: 105°C
    - Ramp rate: 2°C/秒 (緩やか推奨)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ4: First Strand cDNA精製 (オプション)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ※ 通常は精製せず、直接Second strand合成へ進む
  ※ RNA量が非常に多い場合のみ精製推奨

  AMPure XP精製 (実施する場合):
    1. 反応液を室温平衡化: 5分
    2. AMPure XP beads: 40 μL (2.0× ratio)
    3. 混合、インキュベート、洗浄 (前述プロトコル同様)
    4. 溶出: Nuclease-free water 20 μL

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
出力:
  First Strand cDNA: 20 μL
  濃度: 1-5 ng/μL (元のRNA量依存)
  保存: -20°C (短期)、-80°C (長期)
  次ステップ: Second Strand合成へ即座に進む (推奨)
```

**トラブルシューティング**

```yaml
問題1: cDNA収量が低い (<20% conversion)
  原因:
    - RNA分解が進んでいた (RIN <5)
    - RNase汚染
    - 試薬の劣化

  対処:
    - RNA品質再評価、新鮮なサンプル使用
    - RNase除去徹底
    - 試薬の有効期限確認

問題2: リボソーマルRNA混入多量
  原因:
    - SSR (Strand Specificity Reagent) 不足
    - インキュベーション時間不足

  対処:
    - SSR量を1.2倍に増量
    - 42°C インキュベーション時間延長 (20分)

問題3: ストランド特異性喪失
  原因:
    - dUTP混入 (試薬クロスコンタミ)
    - 温度設定ミス

  対処:
    - 新しい試薬使用
    - サーマルサイクラー校正
```

---

## 5. Second strand合成とdsDNA調製

### 5.1 Second Strand合成プロトコル

**dUTP取り込みによるストランドマーキング**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ1: Second Strand反応系調製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  反応系 (80 μL):
    First Strand cDNA: 20 μL (前ステップ全量)
    NEBNext Second Strand Synthesis Buffer (10X): 8 μL
    NEBNext Second Strand Synthesis Enzyme Mix: 4 μL
      (DNA Pol I, RNase H, dUTP mix含有)
    Nuclease-free water: 48 μL
    ────────────────────────────
    総量: 80 μL

  調製手順:
    1. 新しい1.5 mL LoBind tubeに添加:
       a. First Strand cDNA: 20 μL
       b. Nuclease-free water: 48 μL
       c. Second Strand Synthesis Buffer (10X): 8 μL
       d. Second Strand Synthesis Enzyme Mix: 4 μL

    2. ピペッティング混合: 10回
    3. スピンダウン: 3,000 rpm、5秒

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ2: Second Strand合成反応
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  インキュベーション:
    温度: 16°C
    時間: 60分
    機器: サーマルサイクラーまたは冷却ブロック

  反応の進行:
    - RNase H: First strandのRNAテンプレート分解
    - DNA Polymerase I: Second strand合成
    - dUTP: Second strandに取り込まれる (重要!)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ3: dsDNA精製 (AMPure XP)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  目的: 酵素、dNTP、バッファー除去

  手順:
    1. 反応液を室温平衡化: 5分

    2. AMPure XP beads添加: 144 μL (1.8× ratio)
       → 200 bp以上のdsDNA回収

    3. ピペッティング混合: 10回

    4. 室温インキュベート: 5分

    5. マグネティックスタンド: 5分
       溶液が完全に透明になるまで

    6. 上清を完全除去 (廃棄)

    7. ビーズ on stand で80% エタノール洗浄:
       a. 新調製80% エタノール 200 μL添加
       b. 30秒静置
       c. 上清除去
       d. 繰り返し (計2回洗浄)

    8. 残存エタノール完全除去:
       - P10ピペットで残液吸引
       - スタンド上で風乾: 5分
       - 注意: 過乾燥させない (ビーズが割れる)

    9. 溶出:
       - スタンドから外す
       - Nuclease-free water: 53 μL添加
       - ピペッティング混合: 10回
       - 室温インキュベート: 2分
       - マグネティックスタンド: 3分
       - 上清50 μL回収 → 新しい1.5 mL LoBind tube

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
出力:
  dsDNA (dUTP含有Second strand): 50 μL
  濃度: 2-10 ng/μL (元のRNA量依存)
  次ステップ: End Prep/dA-tailing
```

### 5.2 End Prep / dA-tailing

**MinION adapter ligationのための末端調製**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ1: End Prep反応系調製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  使用試薬: NEBNext Ultra II End Repair/dA-Tailing Module
  (または SQK-LSK114キット付属 End Prep試薬)

  反応系 (65 μL):
    dsDNA: 50 μL (前ステップ全量)
    Ultra II End Prep Reaction Buffer: 7 μL
    Ultra II End Prep Enzyme Mix: 3 μL
    ────────────────────────────
    総量: 60 μL

  調製手順:
    1. 0.2 mL PCRチューブに添加:
       a. dsDNA: 50 μL
       b. End Prep Buffer: 7 μL
       c. End Prep Enzyme Mix: 3 μL

    2. ピペッティング混合: 10回
    3. スピンダウン: 3,000 rpm、5秒

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ2: End Prep反応
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  サーマルサイクラー設定:
    ステップ1: 20°C - 5分 (End repair)
    ステップ2: 65°C - 5分 (dA-tailing)
    ステップ3: 4°C - ホールド

  反応内容:
    - 3' overhang除去 (exonuclease活性)
    - 5' overhang充填 (polymerase活性)
    - Blunt end化
    - 3'末端へのdAMP付加 (A-tailing)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ3: AMPure XP精製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  手順:
    1. AMPure XP beads: 60 μL (1.0× ratio)
    2. ピペッティング混合、インキュベート
    3. マグネティックスタンド
    4. 上清除去
    5. 80% エタノール洗浄: 2回 (各200 μL)
    6. 風乾: 5分
    7. 溶出: Nuclease-free water 31 μL
    8. 上清30 μL回収

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
出力:
  End-prepped dsDNA: 30 μL
  濃度: 3-12 ng/μL
  A-overhang: 3'末端
  次ステップ: MinION adapter ligation
```

---

## 6. MinION RNAライブラリ調製 (SQK-RNA002)

### 6.1 Direct RNA Sequencing概要

**参考情報 (cfRNAには非推奨)**

```yaml
SQK-RNA002 Direct RNA Sequencing Kit:

  利点:
    - RNA modification検出 (m6A, m5C, pseudouridine)
    - PCR-freeでバイアス最小
    - ネイティブRNA情報保持

  要求仕様:
    RNA量: 500 ng
    RNA品質: RIN >8
    Poly(A) tail: 必須
    RNA長: >200 nt

  プロトコル概要:
    1. Poly(A) RNA精製 (oligo(dT) beads)
    2. RT adapter ligation (T4 RNA Ligase)
    3. Reverse transcription (RTA付加)
    4. MinION adapter ligation
    5. フローセルローディング

  所要時間: 3-4時間

  PMDA病原体検出への適用評価:
    ✗ 不適
    理由:
      - cfRNA量不足 (10-50倍不足)
      - ウイルスRNAの多くはPoly(A) tailなし
      - スループット低すぎ (0.5-2 Gb)
      - コストパフォーマンス悪い

結論:
  Direct RNA-Seqは研究用途のみ
  PMDA病原体スクリーニングにはcDNA-Seqアプローチ使用
```

---

## 7. cDNA-Seq用ライブラリ調製 (SQK-LSK114)

### 7.1 cDNAライブラリ調製の全体フロー

**cfRNA → MinIONライブラリ完成までの統合プロトコル**

```yaml
全体フロー (1日で完結可能):

  Day 1 - 午前:
    08:00-09:00  cfRNA抽出 (Zymo Kit)
    09:00-10:00  DNase処理
    10:00-11:30  First Strand cDNA合成
    11:30-13:00  Second Strand合成
    13:00-13:30  昼休憩

  Day 1 - 午後:
    13:30-14:30  End Prep / dA-tailing
    14:30-15:00  Adapter ligation (MinION)
    15:00-16:00  AMPure精製、QC
    16:00-16:30  フローセルプライミング
    16:30-17:00  ライブラリローディング、シーケンス開始

  Day 2-3:
    シーケンス実行 (48時間、無人モニタリング)

  Day 4:
    データ解析開始
```

### 7.2 Adapter Ligation (SQK-LSK114)

**cDNAへのMinION adapter付加**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ1: Adapter Ligation反応系調製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  使用試薬: SQK-LSK114 Ligation Sequencing Kit

  反応系 (70 μL):
    End-prepped cDNA: 30 μL (前ステップ全量)
    NEBNext Quick Ligation Reaction Buffer (5X): 25 μL
    Adapter Mix II (AMX II): 5 μL
    Quick T4 DNA Ligase: 10 μL
    ────────────────────────────
    総量: 70 μL

  調製手順:
    1. 新しい1.5 mL DNA LoBind tubeに添加:
       a. End-prepped cDNA: 30 μL
       b. Ligation Buffer (5X): 25 μL
          ※ 粘性高い、P200でゆっくりピペッティング
       c. Adapter Mix II: 5 μL
          ※ 使用直前にボルテックス、スピンダウン
       d. Quick T4 DNA Ligase: 10 μL

    2. ピペッティング混合: 10回 (優しく)
    3. スピンダウン: 3,000 rpm、5秒

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ2: Ligation反応
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  インキュベーション:
    温度: 室温 (20-25°C)
    時間: 10分
    撹拌: 不要 (静置)

  重要:
    - 氷上に置かない (Ligase活性低下)
    - タイマー厳守 (over-ligationでアダプターダイマー増加)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ステップ3: 2段階AMPure XP精製
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  目的:
    - 未反応Adapter除去
    - アダプターダイマー除去
    - 短鎖DNA (<500 bp) 除去

  第1回精製:
    1. AMPure XP beads: 40 μL (0.6× ratio)
       → >700 bp DNA回収、アダプターダイマー除去

    2. ピペッティング混合: 10回
    3. 室温インキュベート: 5分
    4. マグネティックスタンド: 5分
    5. 上清を新しい1.5 mL tubeへ移す (重要!)
       ※ ビーズは捨てない (後で使用)

  第2回精製 (上清を使用):
    6. 上清に追加AMPure XP beads: 40 μL
       → 総beads量: 80 μL = 1.2× ratio (元の70 μL反応系比)

    7. ピペッティング混合: 10回
    8. 室温インキュベート: 5分
    9. マグネティックスタンド: 5分
    10. 上清を完全除去 (廃棄)

    11. SFB (Short Fragment Buffer) 洗浄:
        a. ビーズ on stand、SFB 250 μL添加
        b. 30秒静置
        c. 上清完全除去
        d. 繰り返し (計1回のみ)

    12. LFB (Long Fragment Buffer) 洗浄:
        a. スタンドから外す
        b. LFB 200 μL添加
        c. ピペッティング懸濁: 10回
        d. マグネティックスタンド: 3分
        e. 上清完全除去
        f. 繰り返し (計2回)

    13. 残液完全除去、風乾: 30秒
        ※ 過乾燥厳禁 (最大1分)

    14. 溶出:
        - スタンドから外す
        - Elution Buffer (EB): 15 μL添加
        - ピペッティング混合: 10回
        - 室温インキュベート: 10分
        - マグネティックスタンド: 3分
        - 上清12 μL回収 → 新しい1.5 mL LoBind tube

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
出力:
  MinION-ready cDNAライブラリ: 12 μL
  濃度: 10-30 ng/μL (期待値)
  保存: 4°C (当日ローディング) または -20°C (翌日以降)
```

### 7.3 cDNAライブラリQC

**QC測定項目**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
QC 1: Qubit定量
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  測定:
    Qubit dsDNA HS Assay
    サンプル: ライブラリ 1 μL

  合格基準:
    濃度: 10-50 ng/μL
    総量: 120-600 ng (12 μL全量換算)
    最低許容: 80 ng

  評価:
    優良: >30 ng/μL (>360 ng total)
    良好: 20-30 ng/μL (240-360 ng total)
    許容: 10-20 ng/μL (120-240 ng total)
    不合格: <10 ng/μL (<120 ng total)
      → ライブラリ再調製検討

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
QC 2: TapeStation/Bioanalyzer サイズ分布
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  測定:
    Agilent High Sensitivity D5000 ScreenTape
    または Bioanalyzer High Sensitivity DNA Kit
    サンプル: ライブラリ 1 μL

  合格基準:
    Read N50: >5 kb (cDNAの場合)
    ピークサイズ: 1-10 kb (broad distribution)
    アダプターダイマー (<200 bp): <5%

  評価:
    優良: N50 >8 kb, ダイマー <2%
    良好: N50 5-8 kb, ダイマー 2-5%
    許容: N50 3-5 kb, ダイマー 5-10%
    要改善: N50 <3 kb, ダイマー >10%
      → AMPure精製追加実施

  電気泳動パターン:
    正常:
      ┌────────────────┐
      │         ╱╲     │ ← 1-10 kb broad peak
      │        ╱  ╲    │
      │    ___╱    ╲___│
      │___╱            │
      └────────────────┘
      100bp    5kb   10kb

    異常 (アダプターダイマー多):
      ┌────────────────┐
      │ ╱╲             │ ← 150-200 bp sharp peak (ダイマー)
      │╱  ╲    ╱╲      │ ← ライブラリピーク小さい
      │     ╲__╱  ╲    │
      │             ╲_ │
      └────────────────┘
      100bp    5kb   10kb

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
QC 3: Molarity計算 (フローセルローディング用)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  計算式:
    Molarity (fmol/μL) = (DNA ng/μL × 10^6) / (Avg size bp × 650)

  例: 25 ng/μL, N50 = 6 kb (6,000 bp)
    Molarity = (25 × 10^6) / (6,000 × 650)
             = 25,000,000 / 3,900,000
             = 6.41 fmol/μL

  ローディング量計算:
    目標: 75 fmol (新品フローセル)
    必要体積: 75 / 6.41 = 11.7 μL

  合格基準:
    モル濃度: 3-10 fmol/μL
    必要体積: 8-25 μL (12 μL全量内で足りる)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
QC 4: qPCR品質確認 (オプション)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  目的:
    - 増幅可能なライブラリ確認
    - アダプター正常付加確認

  プライマー:
    MinION adapter特異的プライマー使用
    (Oxford Nanopore提供またはカスタム設計)

  合格基準:
    Ct値: 15-20 (10 ng template使用時)
    融解曲線: 単一ピーク (80-85°C)
```

**QC記録シート**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
cDNAライブラリQC記録
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ライブラリID: _____________________________________
サンプルID: _______________________________________
調製日: YYYY-MM-DD
測定者: _________________ 署名: __________

【Qubit測定】
  測定日時: YYYY-MM-DD HH:MM
  Assay: □ dsDNA HS
  濃度: _______ ng/μL
  総量: _______ ng (12 μL換算)
  判定: □ 優良 □ 良好 □ 許容 □ 不合格

【TapeStation/Bioanalyzer】
  測定日時: YYYY-MM-DD HH:MM
  システム: □ TapeStation D5000  □ Bioanalyzer

  結果:
    Read N50: _______ bp
    ピークサイズ範囲: _______ - _______ bp
    アダプターダイマー: _______ % (<200 bp画分)

  判定: □ 優良 □ 良好 □ 許容 □ 要改善

  電気泳動イメージ添付: □

【モル濃度計算】
  DNA濃度: _______ ng/μL
  平均サイズ: _______ bp
  モル濃度: _______ fmol/μL

  ローディング計算 (75 fmol目標):
    必要体積: _______ μL
    判定: □ 適切 (8-25 μL範囲内)
          □ 要希釈 (>25 μL必要)
          □ 不足 (<8 μL、低精度)

【総合判定】
  □ フローセルローディング可
  □ 条件付き可 (理由: ________________________)
  □ 不可、再調製必要 (理由: __________________)

【次のステップ】
  □ フローセルプライミングへ進む
  □ ライブラリ希釈実施
  □ ライブラリ再調製

確認者: _________________ 署名: __________
確認日時: YYYY-MM-DD HH:MM

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

---

## 8. RNA品質評価とQC基準

### 8.1 RNA抽出後の品質評価

**RIN (RNA Integrity Number) 測定**

```yaml
RIN測定の重要性:
  - RNA分解度の客観的指標
  - RIN 1 (完全分解) 〜 RIN 10 (完全intact)
  - cDNA合成効率に直接影響

測定方法:
  機器: Agilent Bioanalyzer 2100
  キット: RNA 6000 Nano Kit または Pico Kit
  サンプル量: 1 μL (濃度 25-500 ng/μL)

RIN値とcDNA合成への影響:

  RIN 9-10 (Excellent):
    状態: 完全intact、分解なし
    18S/28S rRNA: 明瞭な2本バンド (1:2比)
    cDNA収量: 100% (理論値)
    推奨: そのまま使用

  RIN 7-8 (Good):
    状態: 軽度分解、実用範囲
    18S/28S rRNA: 検出可能だがやや不明瞭
    cDNA収量: 70-90%
    推奨: 使用可能

  RIN 5-6 (Fair):
    状態: 中程度分解
    18S/28S rRNA: 不明瞭、バックグラウンド増加
    cDNA収量: 40-60%
    推奨: 使用可能だが短鎖リード中心
    対策: cDNA合成量増加で補償

  RIN 3-4 (Poor):
    状態: 高度分解
    18S/28S rRNA: 消失、broad smear
    cDNA収量: 20-40%
    推奨: 使用可能だが低品質
    対策: 2倍量のRNA使用

  RIN 1-2 (Degraded):
    状態: ほぼ完全分解
    cDNA収量: <20%
    推奨: 使用非推奨、新しいサンプル推奨

cfRNAの特殊性:
  - cfRNAは本質的に断片化 (自然分解産物)
  - RIN値: 通常2-5 (低くても正常)
  - RINよりもDV200 (200 nt以上の画分%) が重要指標

  DV200基準:
    - >70%: 優良
    - 50-70%: 良好
    - 30-50%: 許容
    - <30%: 要注意
```

### 8.2 cDNA合成後のQC基準

**cDNA conversion効率評価**

```yaml
理論的conversion効率:
  元のRNA 100 ng → 期待cDNA 80-100 ng
  (First + Second strand合成で約2倍、ロス考慮で0.8-1.0×)

実測conversion効率:

  優良 (80-100%):
    - 技術的に理想的
    - 全ステップが最適

  良好 (60-80%):
    - 実用的な範囲
    - 標準的な結果

  許容 (40-60%):
    - やや低いが使用可能
    - RNA品質やや低い可能性

  不良 (<40%):
    - 何らかの問題あり
    - 原因調査と対策必要

計算例:
  元のRNA: 50 ng
  cDNA (Qubit測定): 35 ng
  Conversion効率: 35/50 = 70% (良好)

影響因子:
  - RNA品質 (RIN, DV200)
  - RNase汚染
  - 試薬の鮮度
  - 温度制御の正確性
  - ピペッティング技術
```

### 8.3 シーケンス品質予測

**ライブラリQCからのシーケンス品質予測**

```yaml
予測モデル:

  Input RNA品質 + cDNA conversion + Library QC
    ↓
  予測シーケンス品質

具体例:

  シナリオ1: 最適条件
    RNA: RIN 7, 50 ng
    cDNA conversion: 75% → 37.5 ng
    Library: N50 8 kb, 250 ng total
    予測:
      - Read数: 2-3M reads
      - Total yield: 20-25 Gb
      - Mean Q: Q25-28
      - Duplex rate: 50-60%

  シナリオ2: 標準条件
    RNA: RIN 5, 30 ng
    cDNA conversion: 60% → 18 ng
    Library: N50 5 kb, 150 ng total
    予測:
      - Read数: 1.5-2M reads
      - Total yield: 12-18 Gb
      - Mean Q: Q22-25
      - Duplex rate: 40-50%

  シナリオ3: 低品質条件
    RNA: RIN 3, 20 ng
    cDNA conversion: 40% → 8 ng
    Library: N50 3 kb, 80 ng total
    予測:
      - Read数: 0.8-1.2M reads
      - Total yield: 5-10 Gb
      - Mean Q: Q18-22
      - Duplex rate: 20-30%
      注意: 病原体検出感度低下の可能性
```

---

## 9. RNA特異的トラブルシューティング

### 9.1 RNA分解問題

**症状: RIN <5, DV200 <30%**

```yaml
原因分析:

  原因1: サンプル採取・保存不良
    - 血漿分離遅延 (>1時間)
    - 室温放置
    - 凍結融解繰り返し

  対処:
    - サンプル採取プロトコル再徹底
    - 30分以内の血漿分離
    - -80°C急速凍結
    - アリコート保存 (凍結融解1回のみ)

  原因2: 抽出時のRNase混入
    - ピペット・チューブのRNase汚染
    - 試薬へのRNase混入
    - 氷のRNase汚染

  対処:
    - RNaseZAP処理: 全器具
    - 専用RNase-freeピペット使用
    - 新しい試薬ロット使用
    - RNase-free ice作成 (専用製氷)

  原因3: 抽出プロトコル不適
    - 溶出時の長時間室温放置
    - ビーズ乾燥不足 (エタノール残存)

  対処:
    - 溶出後即座に次ステップへ
    - または -80°C保存
    - ビーズ乾燥時間最適化 (3-5分)
```

### 9.2 cDNA合成低効率問題

**症状: cDNA収量 <40% conversion**

```yaml
原因分析:

  原因1: RNA templateの二次構造
    - GC-richウイルスRNA
    - 強い二次構造形成

  対処:
    - RNA変性処理追加:
      65°C, 5分 → 氷上
    - Primer annealingステップ延長
    - DMSO添加 (最終濃度5%, 二次構造解消)

  原因2: Reverse transcriptase活性不足
    - 試薬劣化
    - 保存温度不適
    - 凍結融解繰り返し

  対処:
    - 新しい試薬ロット使用
    - 酵素アリコート保存 (-80°C)
    - 使用時に初めて解凍

  原因3: RNase汚染
    - cDNA合成中のRNA template分解

  対処:
    - RNaseOUT濃度増加 (2倍量)
    - 反応系のRNase除去徹底

  原因4: Inhibitor混入
    - 血漿由来inhibitor (ヘモグロビンなど)
    - エタノール残存

  対処:
    - RNA精製の徹底 (AMPure XP追加洗浄)
    - 希釈してinhibitor濃度低下
```

### 9.3 ストランド特異性喪失

**症状: (+)鎖と(-)鎖の区別不可**

```yaml
原因分析:

  原因1: Second strand合成時のdUTP混入不足
    - dUTP mix劣化
    - 酵素mix取り違え

  対処:
    - 新しいNEBNext Ultra II Kit使用
    - 試薬ラベル再確認
    - Positive control実験 (既知RNA)

  原因2: UDG処理不適 (もしプロトコルに含む場合)
    - UDG処理時間不足
    - UDG酵素失活

  対処:
    - UDG処理条件最適化
    - 新しいUDG酵素使用

  原因3: PCR増幅時のバイアス (PCR使用時)
    - PCRサイクル数過剰
    - 非特異的増幅

  対処:
    - PCRサイクル数最小化 (<12サイクル)
    - High-fidelity polymerase使用
    - または PCR-free protocol (MinION推奨)

検証方法:
  既知の(+)鎖RNAウイルス (例: PRRSV) をスパイク
  シーケンス後、(+)鎖/(-)鎖リード比を確認
  理想: >95% (+)鎖リード
  問題あり: <80% (+)鎖リード
```

### 9.4 アダプターダイマー過剰

**症状: TapeStationで<200 bpピーク優勢**

```yaml
原因分析:

  原因1: cDNA量不足
    - Adapter過剰 (cDNA:adapter比不適)
    - cDNA合成失敗

  対処:
    - cDNA量を2倍に増やす
    - Adapter量減少 (75% → 50%)
    - cDNA合成プロトコル見直し

  原因2: Ligation時間過剰
    - Over-ligation
    - Self-ligation増加

  対処:
    - Ligation時間短縮 (10分 → 5分)
    - Ligase量減少

  原因3: AMPure XP精製不十分
    - 0.6× ratioでダイマー除去不足

  対処:
    - AMPure比率調整: 0.6× → 0.5×
      (より厳格な >1kb選択)
    - SFB洗浄回数増加 (1回 → 2回)
    - 追加AMPure精製ラウンド実施

  原因4: 高濃度adapter使用
    - Adapter MixⅡ濃度過剰

  対処:
    - Adapter希釈 (50-75%濃度)
```

---

## 10. ALCOA+準拠記録

### 10.1 cDNA合成記録シート

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
cDNA合成記録
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

実施日: YYYY-MM-DD
実施者: [フルネーム] _____________ 署名: _______

【1. 開始RNA情報】
  サンプルID: ___________________________________
  RNA抽出日: YYYY-MM-DD
  RNA濃度 (Qubit): _______ ng/μL
  RNA総量: _______ ng
  RIN値: _______ (Bioanalyzer)
  DV200: _______ %
  RNA品質判定: □ Excellent □ Good □ Fair □ Poor

【2. DNase処理】
  実施時刻: HH:MM
  使用RNA量: _______ μL (_______ ng)

  試薬ロット番号:
    TURBO DNase: _______________
    RNaseOUT: _______________

  インキュベーション:
    開始時刻: HH:MM
    終了時刻: HH:MM
    温度確認: □ 37°C

  精製方法: □ AMPure XP  □ DNase Inactivation Reagent

  結果:
    処理後RNA濃度: _______ ng/μL
    回収率: _______ %
    残存DNA確認: □ qPCR実施 (ΔCt: _______)
                □ 未実施

【3. First Strand cDNA合成】
  開始時刻: HH:MM

  試薬ロット番号:
    NEBNext Ultra II Kit: _______________
    Random Primers: _______________
    Strand Specificity Reagent: _______________

  サーマルサイクラー:
    機器ID: _______________
    プログラム確認: □ 65°C 5min → 25°C 10min →
                      42°C 15min → 70°C 15min

  完了時刻: HH:MM

【4. Second Strand合成】
  開始時刻: HH:MM

  試薬ロット番号:
    Second Strand Synthesis Mix: _______________

  インキュベーション:
    温度: □ 16°C
    時間: □ 60分
    開始時刻: HH:MM
    完了時刻: HH:MM

  精製 (AMPure XP):
    Beads lot: _______________
    溶出体積: _______ μL

【5. End Prep / dA-tailing】
  開始時刻: HH:MM

  試薬ロット番号:
    End Prep Enzyme Mix: _______________

  サーマルサイクラー:
    プログラム確認: □ 20°C 5min → 65°C 5min

  完了時刻: HH:MM

【6. Adapter Ligation】
  開始時刻: HH:MM

  試薬ロット番号:
    Adapter Mix II: _______________
    Quick T4 DNA Ligase: _______________

  Ligation条件:
    温度: □ 室温 (______°C)
    時間: □ 10分

  精製 (2段階AMPure):
    1st round (0.6×): □ 完了
    2nd round (1.2×): □ 完了
    SFB wash: □ 1回
    LFB wash: □ 2回
    溶出体積: _______ μL

  完了時刻: HH:MM

【7. QC結果】
  Qubit測定:
    濃度: _______ ng/μL
    総量: _______ ng

  TapeStation:
    Read N50: _______ bp
    アダプターダイマー: _______ %

  Molarity計算:
    モル濃度: _______ fmol/μL
    75 fmol必要体積: _______ μL

  総合判定:
    □ 優良  □ 良好  □ 許容  □ 不合格

【8. 保存】
  保存場所: ___________________________________
  保存温度: □ 4°C  □ -20°C  □ -80°C
  保存開始日時: YYYY-MM-DD HH:MM

【9. 特記事項】
  ________________________________________________
  ________________________________________________

【10. 確認・承認】
  記録者: [フルネーム] _____________ 署名: _______
  日時: YYYY-MM-DD HH:MM

  確認者: [フルネーム] _____________ 署名: _______
  日時: YYYY-MM-DD HH:MM

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### 10.2 電子データ管理

**cDNA合成データの長期保存**

```yaml
保存対象ファイル:

  1. Bioanalyzer/TapeStation データ:
     - RNA品質 (RIN測定結果)
     - cDNA品質 (N50、サイズ分布)
     - 電気泳動イメージ
     保存期間: 5年以上

  2. Qubit測定データ:
     - RNA濃度
     - cDNA濃度
     - 各精製ステップ濃度
     保存期間: 5年以上

  3. qPCR検証データ (実施時):
     - DNase処理効率
     - ストランド特異性検証
     - Ct値、融解曲線
     保存期間: 5年以上

  4. サーマルサイクラーログ:
     - 温度プロファイル記録
     - 機器ID、日時
     保存期間: 5年以上

ファイル命名規則:
  [SampleID]_[Step]_[Date]_[Replicate].ext

  例:
    PigPlasma001_RNA_QC_20250308_R1.pdf
    PigPlasma001_cDNA_Qubit_20250308.csv
    PigPlasma001_cDNA_TapeStation_20250308.pdf

バックアップ:
  - 本番サーバー (リアルタイム)
  - バックアップNAS (日次)
  - クラウドストレージ (週次)
```

---

## まとめ

本章では、MinION用RNAライブラリ調製の全プロセスを詳述しました。

**重要ポイント:**

1. **cDNA-Seqアプローチの選択**
   - cfRNA量不足のためDirect RNA-Seq不可
   - ストランド特異的cDNA合成必須
   - NEBNext Ultra II Directional Kit推奨

2. **RNase汚染防止の徹底**
   - RNaseZAP処理
   - RNase-free環境
   - 迅速な処理 (RNA抽出→cDNA合成連続実施)

3. **DNase処理の重要性**
   - DNA残存 <1%
   - qPCR検証推奨
   - RNA回収率 >80%

4. **ストランド特異性の確保**
   - dUTP法使用
   - (+)鎖/(-)鎖区別
   - ウイルス複製活性評価可能

5. **品質管理の多層化**
   - RIN/DV200 (RNA品質)
   - Conversion効率 (cDNA合成)
   - N50/アダプターダイマー (ライブラリ品質)

6. **ALCOA+準拠記録**
   - 全ステップ詳細記録
   - 試薬ロット番号トレーサビリティ
   - 電子データ長期保存

**次章予告:**

第7章では、ライブラリ品質評価とDuplex対応について詳述します。Duplex basecallingの原理、最適化戦略、Q30+データ取得方法について解説します。

---

**文書情報**
- 作成日: 2025-03-08
- バージョン: 1.0
- 作成者: MinIONメタゲノム解析プロトコル開発チーム
- 承認者: [承認者名]
- 次回改訂予定: 2025-09-08
