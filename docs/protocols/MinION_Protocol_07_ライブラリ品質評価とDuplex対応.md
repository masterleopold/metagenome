# MinION用メタゲノム解析プロトコル
# 第7章: ライブラリ品質評価とDuplex対応

## 目次

1. [ライブラリ品質評価の重要性](#1-ライブラリ品質評価の重要性)
2. [Duplex Sequencingの原理と利点](#2-duplex-sequencingの原理と利点)
3. [Duplex対応ライブラリ調製最適化](#3-duplex対応ライブラリ調製最適化)
4. [ライブラリQC詳細プロトコル](#4-ライブラリqc詳細プロトコル)
5. [Duplex basecalling設定と最適化](#5-duplex-basecalling設定と最適化)
6. [Q30+データ取得戦略](#6-q30データ取得戦略)
7. [ライブラリ保存と安定性](#7-ライブラリ保存と安定性)
8. [バッチ間変動管理](#8-バッチ間変動管理)
9. [ライブラリQCトラブルシューティング](#9-ライブラリqcトラブルシューティング)
10. [ALCOA+準拠記録](#10-alcoa準拠記録)

----

## 1. ライブラリ品質評価の重要性

### 1.1 ライブラリ品質がシーケンスに与える影響

**品質-スループット-精度の三角関係**

```yaml
高品質ライブラリの特徴:
  1. 適切なサイズ分布:
     - 長鎖DNA優位 (N50 >10 kb)
     - 短鎖DNA minimal (<500 bp <5%)
     - 影響: Read長、カバレッジ均一性

  2. 高純度:
     - Protein汚染なし (A260/A280 >1.8)
     - 塩類minimal
     - Inhibitor除去完全
     - 影響: ポア寿命、シーケンス速度

  3. 適切なアダプター付加:
     - Ligation効率 >80%
     - アダプターダイマー <5%
     - 影響: アクティブチャンネル数、データ量

  4. DNA完全性:
     - ニックなし
     - 一本鎖DNA minimal
     - 影響: Duplex率、精度

品質とシーケンス結果の相関:

  優良ライブラリ (N50 >12 kb, ダイマー <2%):
    → Total yield: 25-30 Gb
    → Duplex rate: 50-60%
    → Mean Q score: Q28-30
    → Active channel持続: >36時間

  標準ライブラリ (N50 8-12 kb, ダイマー 2-5%):
    → Total yield: 18-25 Gb
    → Duplex rate: 40-50%
    → Mean Q score: Q25-28
    → Active channel持続: 24-36時間

  低品質ライブラリ (N50 <8 kb, ダイマー >5%):
    → Total yield: <15 Gb
    → Duplex rate: <30%
    → Mean Q score: Q20-25
    → Active channel持続: <24時間
    → 警告: PMDA病原体検出困難の可能性
```

### 1.2 PMDA病原体検出に必要なライブラリ品質基準

**最低要求仕様**

```yaml
DNA/cDNAライブラリ共通基準:

  必須要件 (Minimum Requirements):
    総量: ≥100 ng
    濃度: ≥10 ng/μL
    Read N50: ≥8 kb
    アダプターダイマー: <10%
    純度 (A260/A280): 1.8-2.0

  推奨基準 (Recommended):
    総量: ≥200 ng
    濃度: ≥20 ng/μL
    Read N50: ≥12 kb
    アダプターダイマー: <5%
    純度 (A260/A280): 1.9-2.0

  最適基準 (Optimal):
    総量: ≥300 ng
    濃度: ≥30 ng/μL
    Read N50: ≥15 kb
    アダプターダイマー: <2%
    純度 (A260/A280): 1.95-2.05

  Duplex対応追加基準:
    DNA integrity: >90% (double-stranded)
    Nick frequency: <1 per 20 kb
    GC content: 検証 (40-60%推奨)

理由:
  91病原体の多様性:
    - GC含量: 25% (AT-rich virus) ～ 70% (GC-rich bacteria)
    - ゲノムサイズ: 1 kb (small virus) ～ 5 Mb (bacteria)
    - 存在量: 0.001-5% (total cfDNA/cfRNA)

  → 高品質ライブラリで低存在量病原体も検出
```

----

## 2. Duplex Sequencingの原理と利点

### 2.1 Duplex Sequencingとは

**DNA二本鎖の両方をシーケンスして精度向上**

```yaml
原理:

  通常のNanoporeシーケンス (Simplex):
    DNA二本鎖 ═══════════════
                  ↓
    一本鎖のみ通過 ───────────── (Template strand)
                  ↓
    Basecall精度: Q20-25 (99-99.9% accuracy)

  Duplexシーケンス:
    DNA二本鎖 ═══════════════
                  ↓
    Hairpin adapterでループ化
                  ↓
    両鎖連続通過 ───────┐
                      ├──── (Template strand)
                      └──── (Complement strand)
                  ↓
    Consensus basecall (両鎖情報統合)
                  ↓
    Basecall精度: Q30+ (99.9%+ accuracy)

Duplex化の仕組み:
  1. Adapter ligation:
     5'─────────────3'  Template strand
     3'─────────────5'  Complement strand
          ↓ Hairpin adapter付加
     5'─────────────┐
                    │ Hairpin loop
     3'─────────────┘

  2. Nanopore通過:
     Template strand → Hairpin → Complement strand
     連続で読まれる

  3. Basecalling:
     両鎖の情報を比較、consensus決定
     → 高精度配列 (Q30+)
```

### 2.2 Duplexの利点とPMDA病原体検出への応用

**精度向上による病原体同定信頼性の劇的改善**

```yaml
Simplex vs Duplex比較:

  パラメータ           Simplex         Duplex
  ────────────────────────────────────────────
  精度 (Q score)      Q20-25          Q30+
  塩基レベル正確性      99-99.9%        99.9%+
  SNP検出信頼性        中              高
  系統解析精度          低-中           高
  偽陽性率             1-5%            <0.1%
  データ量 (同一ラン)   100%            40-60%
  解析時間             1×              1.5-2×

PMDA病原体検出における利点:

  1. 偽陽性の劇的削減:
     Simplex: 1,000 readsで10-50個のエラー
             → 病原体誤同定リスク
     Duplex:  1,000 readsで1-5個のエラー
             → 高信頼性同定

  2. 低存在量病原体の確実な検出:
     - 0.01% (1/10,000) の病原体DNA
     - Simplexでは偽陽性と区別困難
     - Duplexで確実に識別可能

  3. 薬剤耐性遺伝子の検出:
     - 1塩基変異の正確検出
     - 耐性株の確実な同定

  4. ウイルス系統解析:
     - SNP levelの系統分類
     - ワクチン株 vs 野生株識別

  5. 複数病原体の同時検出:
     - 混合感染の正確な判定
     - 各病原体の正確な定量

推奨:
  PMDA 91病原体スクリーニングでは
  Duplex basecalling必須
```

### 2.3 Duplex率を決定する因子

**ライブラリ品質とシーケンス条件の影響**

```yaml
Duplex率 = (Duplex reads / Total reads) × 100

影響因子:

  1. DNA/cDNA完全性 (最重要):
     DNA integrity >95%: Duplex rate 50-60%
     DNA integrity 80-95%: Duplex rate 30-50%
     DNA integrity <80%: Duplex rate <30%

  2. DNA長:
     Average length >10 kb: Duplex rate 50-60%
     Average length 5-10 kb: Duplex rate 40-50%
     Average length <5 kb: Duplex rate <30%

  3. Adapter ligation効率:
     Ligation >90%: Duplex rate 50-60%
     Ligation 70-90%: Duplex rate 35-50%
     Ligation <70%: Duplex rate <35%

  4. シーケンス深度:
     >2×カバレッジ: Duplex rate 50-60%
     1-2×カバレッジ: Duplex rate 30-40%
     <1×カバレッジ: Duplex rate <20%
     理由: 同じ分子が複数回読まれる必要

  5. フローセル品質:
     Active pores >1,500: Duplex rate 50-60%
     Active pores 1,000-1,500: Duplex rate 40-50%
     Active pores <1,000: Duplex rate <40%

  6. ラン時間:
     48-72時間: Duplex rate 50-60%
     24-48時間: Duplex rate 40-50%
     <24時間: Duplex rate <35%

最適化戦略:
  - DNA抽出時の物理的シェアリング最小化
  - AMPure精製でのピペッティング優しく
  - ライブラリ保存温度管理 (-20°C)
  - 48時間シーケンスラン
  - 高品質フローセル使用 (>1,500 pores)
```

----

## 3. Duplex対応ライブラリ調製最適化

### 3.1 DNA抽出時の最適化

**長鎖・高完全性DNA回収**

```yaml
cfDNA/cfRNA抽出最適化:

  目標:
    - 断片化最小化
    - 長鎖DNA保持
    - 一本鎖DNA minimal

  Zymo Quick-cfDNA Kit使用時の最適化:

    1. サンプル取り扱い:
       - 血漿の凍結融解: 1回のみ (複数回厳禁)
       - 解凍方法: 氷上でゆっくり (室温急速解凍禁止)
       - ピペッティング: 非常に優しく (10回以下)

    2. 溶出条件最適化:
       - 溶出バッファー温度: 55°C (通常50°C)
       - インキュベーション時間: 5分 → 10分
       - 理由: 長鎖DNAの完全溶出

    3. 磁気ビーズ処理:
       - ボルテックス厳禁 (タッピングまたは反転)
       - マグネティックスタンド時間: 十分に (5-10分)
       - エタノール洗浄: P1000で優しく

    4. 濃縮方法 (低濃度時):
       - 真空濃縮器使用 (SpeedVac)
       - 設定: 30°C、15分
       - 過乾燥注意 (DNA変性リスク)
```

### 3.2 AMPure XP精製の最適化

**サイズ選択とDNA完全性の両立**

```yaml
AMPure XP比率の最適化:

  目的別推奨比率:

    長鎖DNA最大化 (Duplex率優先):
      1st cleanup: 0.6× (>1 kb選択)
      2nd cleanup: なし
      メリット: 長鎖DNA最大保持
      デメリット: アダプターダイマー残存可能性

    バランス型 (推奨):
      1st cleanup: 0.5× (>2 kb選択)
      2nd cleanup: 0.4× 追加 (>3 kb選択)
      メリット: 長鎖保持 + ダイマー除去
      デメリット: 若干のロス

    高純度型:
      1st cleanup: 0.4× (>3 kb選択)
      2nd cleanup: 0.6× (短鎖完全除去)
      メリット: アダプターダイマー完全除去
      デメリット: 長鎖もやや除去

  PMDA病原体検出推奨:
    バランス型 (0.5× + 0.4×)
    理由: 長鎖保持 + 低バックグラウンド

ビーズ洗浄最適化:

  エタノール洗浄:
    濃度: 80% (freshly prepared)
    回数: 2回
    各洗浄時間: 30秒
    注意: 過剰洗浄でDNA変性

  乾燥時間:
    室温: 3-5分
    警告サイン (過乾燥):
      - ビーズ表面のひび割れ
      - 茶色化
      → DNA回収率激減

    対処:
      - 目視で湿り気残存確認
      - 過乾燥した場合: 溶出時間2倍

  溶出条件:
    バッファー: Nuclease-free water または EB
    温度: 室温 (50°Cは変性リスク)
    インキュベーション: 5-10分
    混合: ピペッティング10回 (優しく)
```

### 3.3 Adapter Ligation最適化

**高効率・均一なアダプター付加**

```yaml
Ligation効率向上戦略:

  1. DNA末端調製の徹底:
     End Prep反応:
       - 温度精度: ±0.5°C以内
       - 時間厳守: 5分ジャスト (過剰反応でDNA損傷)
       - A-tailing効率: >95%

     確認方法:
       ゲル電気泳動で+25 bp シフト確認
       (dAMP付加による)

  2. Adapter:DNA ratio最適化:
     理想比: Adapter 10-20 fmol : DNA 1 μg
     過剰Adapter: ダイマー増加
     不足Adapter: Ligation効率低下

     計算例:
       DNA: 200 ng, 10 kb平均
       DNA molarity: (200×10^6)/(10,000×650) = 30.8 fmol
       Adapter: 30.8 × 15 = 462 fmol必要
       Adapter Mix II: 通常5 μLで充分

  3. Ligase活性最適化:
     酵素: NEBNext Quick T4 DNA Ligase
     温度: 20-25°C (室温)
     時間: 10分 (厳守)
     ATP濃度: Ligation Buffer内で最適化済み

     注意:
       - 氷上禁止 (Ligase失活)
       - >15分でover-ligation (ダイマー増)

  4. Ligation Buffer fresh使用:
     - 開封後3ヶ月以内
     - -20°C保管
     - 凍結融解3回まで
     - PEG含有 (粘性高い、ピペッティング注意)

  5. Nick sealing (オプション):
     一本鎖nick残存によるDuplex失敗防止

     追加ステップ:
       Ligation後にNick repair:
         - NEBNext FFPE DNA Repair Mix使用
         - 20°C, 15分
         - AMPure精製

     効果:
       Duplex rate: 40% → 55%に向上
```

----

## 4. ライブラリQC詳細プロトコル

### 4.1 Qubit定量プロトコル

**正確な濃度測定**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Qubit dsDNA HS Assay プロトコル
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

機器: Qubit 4 Fluorometer

試薬:
  - Qubit dsDNA HS Assay Kit (Invitrogen, Q32851)
  - Qubit Assay Tubes (500 μL clear tubes)

手順:

  ステップ1: Working solution調製
    計算: (サンプル数 + 2 standards + 1 extra) × 200 μL

    例: 3サンプル測定
      必要量: (3 + 2 + 1) × 200 = 1,200 μL

    調製:
      Qubit dsDNA HS Buffer: 1,200 μL
      Qubit dsDNA HS Reagent: 6 μL (1:200希釈)
      混合: ボルテックス 5秒

  ステップ2: Standard調製
    Standard #1 (Low):
      Working solution: 190 μL
      Standard #1: 10 μL
      → 最終0 ng/mL

    Standard #2 (High):
      Working solution: 190 μL
      Standard #2: 10 μL
      → 最終10 ng/mL (HS assay)

  ステップ3: Sample調製
    各サンプル:
      Working solution: 199 μL
      Sample: 1 μL
      → 1:200希釈

    注意:
      - ライブラリ濃度が>100 ng/μLの場合
      - さらに希釈必要 (1 μL sample → 0.5 μL + 0.5 μL water)

  ステップ4: インキュベーション
    室温、2分間
    遮光 (アルミホイルで覆う)

  ステップ5: 測定
    Qubit設定:
      - Assay選択: dsDNA HS
      - Sample volume: 1 μL
      - Dilution factor: 200

    測定順序:
      1. Standard #1
      2. Standard #2
      3. Sample 1
      4. Sample 2
      5. Sample 3...

  ステップ6: 結果記録
    測定値: Qubit表示濃度 (自動計算済み)
    単位: ng/μL

    総量計算:
      Total DNA (ng) = 濃度 (ng/μL) × Volume (μL)

      例: 25 ng/μL, 12 μL
        Total = 25 × 12 = 300 ng

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
QC基準:

  優良: 30-50 ng/μL (360-600 ng total)
  良好: 20-30 ng/μL (240-360 ng total)
  許容: 10-20 ng/μL (120-240 ng total)
  不合格: <10 ng/μL (<120 ng total)

注意事項:
  - 測定は3連で実施 (CV <5%)
  - 標準偏差大きい場合: ライブラリ不均一の可能性
  - Qubitは二本鎖DNA特異的 (一本鎖は測定されない)
```

### 4.2 TapeStation/Bioanalyzer プロトコル

**サイズ分布と品質評価**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Agilent TapeStation 4200 - D5000 ScreenTape
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

試薬:
  - D5000 ScreenTape (5067-5588)
  - D5000 Reagents (5067-5589)
  - D5000 Ladder (5067-5590)

サンプル準備:

  ステップ1: 濃度調整
    最適濃度範囲: 5-50 ng/μL
    サンプル量: 2 μL必要

    例: ライブラリ濃度30 ng/μL
      → そのまま使用可能

    例: ライブラリ濃度100 ng/μL
      → 希釈必要
      希釈: 1 μL library + 2 μL water = 33.3 ng/μL

  ステップ2: Sample buffer混合
    Strip tubeに添加:
      D5000 Sample Buffer: 10 μL
      Sample (DNA): 1 μL
      ────────────────────────
      Total: 11 μL

    Ladder準備:
      D5000 Sample Buffer: 10 μL
      D5000 Ladder: 1 μL

  ステップ3: ボルテックスとスピンダウン
    ボルテックス: 2,000 rpm, 1分
    スピンダウン: 1,000×g, 5秒

  ステップ4: TapeStation測定
    1. ScreenTapeとLoading tipsをセット
    2. Strip tubeをセット
    3. ソフトウェアでサンプル情報入力
    4. "Start" をクリック
    5. 測定時間: 約10分 (12サンプル)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
結果解析:

  1. サイズ分布グラフ:
     X軸: DNA size (bp)
     Y軸: Fluorescence intensity (FU)

     正常パターン:
       ┌──────────────────────┐
       │            ╱╲        │
       │          ╱    ╲      │
       │        ╱        ╲    │
       │    ___╱          ╲___│
       │___╱                  │
       └──────────────────────┘
       100  1k  5k  10k  20k (bp)

     異常パターン (アダプターダイマー):
       ┌──────────────────────┐
       │╱╲                    │ ← 150-200 bp sharp peak
       │  ╲     ╱╲            │
       │   ╲___╱  ╲___        │
       │                 ╲___ │
       └──────────────────────┘
       100  1k  5k  10k  20k (bp)

  2. 重要パラメータ:

     Region Summary:
       Size [bp]: ピークサイズ
       % of Total: 各サイズ画分の割合
       Molarity [pmol/l]: モル濃度

     Sample Summary:
       Concentration: 濃度 (ng/μL)
       DIN (DNA Integrity Number): 1-10
         - DIN >7: 優良
         - DIN 5-7: 良好
         - DIN <5: 断片化
       Average Size: 平均サイズ

     N50計算:
       - ソフトウェア自動計算
       - または手動計算:
         50%の累積データ量に達するサイズ

  3. アダプターダイマー定量:
     Region設定: 100-300 bp
     % of Total: 5%未満が目標

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
QC判定基準:

  優良ライブラリ:
    N50: >12 kb
    DIN: >7
    アダプターダイマー: <2%
    ピーク形状: Broad, smooth

  良好ライブラリ:
    N50: 8-12 kb
    DIN: 5-7
    アダプターダイマー: 2-5%
    ピーク形状: やや狭いがacceptable

  許容ライブラリ:
    N50: 5-8 kb
    DIN: 3-5
    アダプターダイマー: 5-10%
    ピーク形状: 短鎖寄り

  不合格:
    N50: <5 kb
    DIN: <3
    アダプターダイマー: >10%
    対応: AMPure追加精製または再調製
```

### 4.3 NanoDrop分光光度計測定

**純度評価**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
NanoDrop One/2000 純度測定
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

測定項目:
  - A260/A280比 (タンパク質汚染評価)
  - A260/A230比 (塩類・有機溶媒汚染評価)
  - 吸光スペクトル (汚染物質検出)

プロトコル:

  ステップ1: Blank測定
    Nuclease-free water: 1.5 μL
    測定台にロード
    "Blank" をクリック

  ステップ2: Sample測定
    ブランク拭き取り (Kimwipe)
    Sample: 1.5 μL
    測定台にロード
    測定タイプ選択: "dsDNA"
    "Measure" をクリック

  ステップ3: 結果記録
    Concentration (ng/μL): DNA濃度
    A260/A280: 純度指標1
    A260/A230: 純度指標2

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
判定基準:

  A260/A280比:
    1.8-2.0: 高純度 (タンパク質汚染minimal)
    1.6-1.8: やや低い (軽度タンパク質汚染)
    <1.6: タンパク質汚染あり
    >2.0: RNA汚染またはアルカリ変性

    対処 (1.6-1.8):
      - AMPure XP追加精製
      - Proteinase K処理

  A260/A230比:
    2.0-2.2: 高純度
    1.5-2.0: やや低い (軽度汚染)
    <1.5: 塩類・EDTA・フェノール汚染

    対処 (<1.5):
      - エタノール沈殿
      - AMPure XP追加精製 (複数回)

  吸光スペクトル評価:
    正常:
      - 260 nmにピーク
      - 230 nm、280 nmは低い
      - Smooth curve

    異常:
      - 230 nmピーク高い: 塩類汚染
      - 270 nm肩あり: フェノール汚染
      - 320-400 nm上昇: 濁り・微粒子

注意:
  NanoDropは総核酸濃度測定 (一本鎖+二本鎖)
  Qubitとの乖離大きい場合:
    - 一本鎖DNA多い
    - RNA混入
    - タンパク質干渉
```

### 4.4 Fragment Analyzer (Advanced QC)

**高精度サイズ分布解析**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Agilent Fragment Analyzer (オプション)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

利点:
  - TapeStationより高分解能
  - 定量精度高い
  - >50 kbのDNA測定可能
  - マルチプレックス対応

使用ケース:
  - 超長鎖DNA解析 (N50 >20 kb)
  - 正確なアダプターダイマー定量
  - 複数サンプル同時QC

キット:
  - DNF-464 HS Large Fragment Kit (1 bp - 50 kb)

測定パラメータ:
  - Size distribution (1 bp resolution)
  - Concentration (高精度)
  - N50, N90, N10
  - GC content
  - Smear analysis

PMDA病原体検出での使用:
  推奨: TapeStationで充分
  Fragment Analyzer: 研究グレード解析時のみ
```

----

## 5. Duplex basecalling設定と最適化

### 5.1 MinKNOW Duplex設定

**Duplex basecallingの有効化**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MinKNOW Duplex Basecalling設定
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

シーケンス開始前設定:

  ステップ1: Basecalling mode選択
    位置: "Output" タブ

    選択肢:
      ○ Super accuracy (Duplex) ← 選択
        - Duplex basecalling有効
        - GPU必須
        - 最高精度 (Q30+)

      ○ High accuracy
        - Simplex basecallのみ
        - 精度: Q20-25
        - PMDA用途には不適

      ○ Fast
        - 低精度、リアルタイム困難
        - 使用非推奨

  ステップ2: Basecaller configuration
    Model自動選択:
      - R9.4.1: dna_r9.4.1_450bps_sup.cfg
      - R10.4.1: dna_r10.4.1_e8.2_400bps_sup.cfg

    Duplex-specific settings:
      ☑ Enable duplex basecalling
      ☑ Generate duplex reads
      Duplex mode: "Automatic"

  ステップ3: Output filtering (オプション)
    Duplex優先出力:
      ☑ Output only duplex reads
        → Simplexリードは出力されない
        → データ量削減、精度最大化

      または
      ☑ Output both simplex and duplex
        → 全リード出力
        → 後処理でDuplex抽出可能

    PMDA推奨: Both出力
      理由: Simplexでも有用情報あり

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GPU設定:

  /opt/ont/minknow/conf/app_conf/guppy_config.toml

  重要パラメータ:
    gpu_runners_per_device = 6
      - Duplex処理に最適な並列数
      - GPU VRAM 24 GBの場合: 6-8推奨
      - GPU VRAM 16 GBの場合: 4-6推奨

    chunks_per_runner = 768
      - Duplex処理用に増量
      - デフォルト512から増量

    chunk_size = 2000
      - シグナル処理単位
      - Duplex最適値: 2000-3000

  適用:
    sudo systemctl restart minknow
```

### 5.2 Dorado Duplex Basecalling (オフライン)

**リアルタイム処理失敗時の対処**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Dorado Standalone Duplex Basecalling
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

使用ケース:
  - MinKNOW basecallingが追いつかない
  - GPU性能不足
  - ラン後に高精度basecalling実施

インストール:
  # Dorado v0.5.0以降
  wget https://cdn.oxfordnanoporetech.com/software/analysis/\
    dorado-0.5.0-linux-x64.tar.gz
  tar -xzf dorado-0.5.0-linux-x64.tar.gz
  export PATH=/path/to/dorado/bin:$PATH

Duplex basecalling実行:

  基本コマンド:
    dorado duplex \
      [model] \
      [pod5_dir] \
      --emit-fastq \
      > output_duplex.fastq

  詳細例:
    dorado duplex \
      dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
      /var/lib/minknow/data/run_20250308/pod5_pass/ \
      --emit-fastq \
      --min-qscore 10 \
      --kit-name SQK-LSK114 \
      --device cuda:0 \
      --verbose \
      > metagenome_duplex_20250308.fastq

  パラメータ説明:
    model: dna_r10.4.1_e8.2_400bps_sup@v4.2.0
      - R10.4.1フローセル用
      - Super accuracy model
      - v4.2.0: 最新モデル

    pod5_dir: POD5ファイルディレクトリ
      - MinKNOWがFAST5→POD5変換済み

    --emit-fastq: FASTQ形式出力
      - または --emit-sam (BAM用)

    --min-qscore 10: 最低Q score filter
      - Duplex: Q30以上が多数
      - Simplex: Q20前後

    --kit-name: ライブラリキット指定
      - SQK-LSK114
      - SQK-RNA002など

    --device cuda:0: GPU指定
      - cuda:0: 1枚目のGPU
      - cuda:all: 全GPU使用

    --verbose: 詳細ログ出力

  出力:
    metagenome_duplex_20250308.fastq
      - Duplexリード + Simplexリード混在
      - Duplexリードは高Q score (Q30+)

  Duplexリードのみ抽出:
    cat metagenome_duplex_20250308.fastq | \
      NanoFilt -q 25 \
      > metagenome_duplex_Q25plus.fastq

    理由: Q25以上のリードは実質Duplex

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
処理時間とリソース:

  データ量: 20 Gb (raw signal)
  GPU: NVIDIA RTX 4090 (24 GB)
  処理時間: 6-10時間 (Duplex basecalling)

  比較:
    Simplex basecalling: 2-3時間
    Duplex basecalling: 6-10時間 (3-5倍時間)

  並列化:
    複数GPUで処理高速化可能
    --device cuda:all
```

----

## 6. Q30+データ取得戦略

### 6.1 Q30+データの重要性

**病原体検出精度への影響**

```yaml
Q scoreと精度の関係:

  Q score定義:
    Q = -10 × log10(P)
    P = error probability

  具体例:
    Q10: 90% accuracy (10% error)
    Q20: 99% accuracy (1% error)
    Q25: 99.68% accuracy (0.32% error)
    Q30: 99.9% accuracy (0.1% error)
    Q40: 99.99% accuracy (0.01% error)

  1,000 bpリードでの期待エラー数:
    Q10: 100 errors
    Q20: 10 errors
    Q25: 3.2 errors
    Q30: 1 error
    Q40: 0.1 errors

PMDA病原体検出でのQ30+必要性:

  シナリオ: 10 kb病原体ゲノム検出

  Q20データ (Simplex):
    10,000 bp × 0.01 error rate = 100 errors/genome
    → SNP判定困難
    → 系統解析不可能
    → 偽陽性リスク高

  Q30データ (Duplex):
    10,000 bp × 0.001 error rate = 10 errors/genome
    → SNP判定可能
    → 系統解析可能
    → 偽陽性リスク低

  推奨:
    総データ量: 20 Gb
    うちQ30+: 8 Gb以上 (40%以上)
    → Duplex rate 40-50%達成必要
```

### 6.2 Q30+データ最大化戦略

**ライブラリからシーケンスまでの統合最適化**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Phase 1: ライブラリ品質最適化
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  目標: Duplex率50%以上

  1. DNA完全性最大化:
     - 抽出時物理的シェアリング回避
     - ピペッティング最小限 (<10回)
     - ボルテックス厳禁
     - 氷上保管徹底

  2. 長鎖DNA優先回収:
     - AMPure比率: 0.4-0.5×
     - 目標N50: >12 kb
     - >20 kb画分: >20%

  3. Nick除去:
     - FFPE DNA Repair Mix使用
     - Nick sealing効率: >95%
     - Duplex率向上: +10-15%

  4. 高効率Adapter ligation:
     - Ligation効率: >90%
     - Adapter:DNA比最適化
     - 時間厳守 (10分)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Phase 2: シーケンス条件最適化
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. 高品質フローセル使用:
     - Active pores: >1,500
     - 新品フローセル推奨
     - 保管期限内 (3ヶ月以内)

  2. 適切なローディング量:
     - 75-100 fmol (1,500+ pores)
     - 過剰ローディング避ける (ポア飽和防止)

  3. 長時間ラン:
     - 48-72時間
     - Duplex生成に時間必要
     - 24時間以降Duplex率上昇

  4. 温度安定性:
     - フローセル温度: 34±0.5°C
     - 室温: 20-25°C
     - 温度変動でポア劣化

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Phase 3: Basecalling最適化
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. Super accuracy model使用:
     - 最新モデル (v4.2.0以降)
     - Duplex mode有効

  2. GPU性能確保:
     - RTX 4090 or A6000推奨
     - リアルタイム処理可能

  3. Duplex calling parameters:
     - Duplex window: デフォルト
     - Pair finding: 自動最適化

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
期待される結果:

  最適条件下:
    Total yield: 25-30 Gb
    Duplex rate: 50-60%
    Q30+ yield: 12-18 Gb
    Mean Q score: Q28-30

  標準条件:
    Total yield: 18-25 Gb
    Duplex rate: 40-50%
    Q30+ yield: 8-12 Gb
    Mean Q score: Q25-28

  両条件ともPMDA病原体検出に充分
```

----

## 7. ライブラリ保存と安定性

### 7.1 ライブラリ保存条件

**長期安定性確保**

```yaml
保存温度と期間:

  短期保存 (当日〜3日):
    温度: 4°C
    容器: 1.5 mL DNA LoBind tube
    光: 遮光 (アルミホイル)
    安定性: >95%活性維持

  中期保存 (1週間〜1ヶ月):
    温度: -20°C
    容器: 1.5 mL DNA LoBind tube
    アリコート: 12 μL × 1本 (凍結融解避ける)
    安定性: >90%活性維持

  長期保存 (1ヶ月〜6ヶ月):
    温度: -80°C
    容器: Cryovial (1.5 mL screw cap)
    アリコート: 6 μL × 2本 (backup含む)
    安定性: >85%活性維持

  超長期保存 (>6ヶ月):
    推奨せず
    理由:
      - Adapter degradation
      - DNA nick増加
      - Duplex率低下

保存時の注意:

  1. 凍結融解回数制限:
     最大2回まで
     3回目以降: 活性50%以下

  2. 解凍方法:
     - 氷上で緩やかに (15-30分)
     - 室温急速解凍禁止
     - 37°C加温厳禁

  3. ラベル記載事項:
     - Sample ID
     - Library prep date
     - Concentration (ng/μL)
     - Volume (μL)
     - Kit used (SQK-LSK114)
     - Freeze-thaw count

  4. 在庫管理:
     - ライブラリ在庫台帳記録
     - 先入先出 (FIFO) 管理
     - 有効期限監視
```

### 7.2 保存安定性の評価

**保存後のQC再実施**

```yaml
保存後QCプロトコル:

  1週間保存後:
    Qubit測定: 濃度確認 (±10%許容)
    TapeStation: N50確認 (±15%許容)
    判定:
      濃度低下<10%: 問題なし
      濃度低下10-20%: 使用可能
      濃度低下>20%: 原因調査

  1ヶ月保存後:
    Qubit測定: 必須
    TapeStation: 必須
    qPCR (オプション): Ligation効率確認
    判定:
      濃度低下<15%: 問題なし
      N50低下<20%: 許容範囲
      Adapter付加率>70%: 使用可能

  3ヶ月以上保存:
    推奨: 新鮮ライブラリ再調製
    理由: Duplex率著明低下のリスク

保存劣化の兆候:

  - 濃度低下>20%
  - N50低下>30%
  - TapeStationで低分子量smear増加
  - シーケンス時:
    * Active channel急速減少
    * Duplex率<20%
    * Mean Q score低下

対処:
  新しいサンプルから再調製
```

----

## 8. バッチ間変動管理

### 8.1 ライブラリ調製のバッチ管理

**再現性確保**

```yaml
バッチ定義:
  1バッチ = 同日・同一ロット試薬で調製したライブラリ群
  推奨バッチサイズ: 1-24サンプル

バッチ管理の重要性:
  - 試薬ロット差の影響排除
  - 操作者技術のばらつき検出
  - システマティックエラー検出
  - ALCOA+準拠 (Consistency)

バッチQC基準:

  Positive Control (必須):
    各バッチに既知DNA (λ DNAなど) 含める
    目標値:
      - 収量: 200±30 ng
      - N50: 20±3 kb
      - Adapter ligation: >85%

  Negative Control (推奨):
    No template control
    目標: 検出限界以下

  バッチ内変動:
    サンプル間CV <20%
    許容範囲: CV 20-30%
    要調査: CV >30%

  バッチ間変動:
    バッチ平均値のCV <15%
    許容範囲: CV 15-25%
    要調査: CV >25%

変動要因分析:

  大きなバッチ間変動の原因:
    - 試薬ロット差
    - 操作者技術差
    - 機器校正ずれ (サーマルサイクラー)
    - 環境条件 (温度、湿度)

  対処:
    - 試薬ロット固定
    - SOP厳守、トレーニング強化
    - 機器定期校正
    - 環境モニタリング
```

### 8.2 内部標準の使用

**定量精度向上**

```yaml
Internal Standard Spike-in:

  目的:
    - 絶対定量の基準
    - ライブラリ調製効率評価
    - バッチ間正規化

  推奨標準品:
    - ERCC RNA Spike-In Mix (RNAライブラリ用)
    - Lambda DNA (DNAライブラリ用)
    - Synthetic pathogen genome (カスタム)

  使用方法:
    1. 既知量のstandard追加:
       例: λ DNA 1 ng/μLを10 μL
       サンプルDNA 100 ngに混合

    2. ライブラリ調製実施
       (standardもライブラリ化される)

    3. シーケンス後、standard mapping:
       λ genome (48.5 kb) への mapping read数

    4. 定量計算:
       Expected reads = (Input λ DNA / Total DNA) × Total reads
       Observed reads = Mapping read count
       Recovery rate = (Observed / Expected) × 100%

    5. サンプル定量補正:
       Corrected量 = 測定量 × (100% / Recovery rate)

  効果:
    - バッチ間補正可能
    - 真の病原体量推定
    - QC fail早期検出

PMDA病原体検出での使用:
  推奨: 研究段階では有用
  臨床検査: validation後に導入検討
```

----

## 9. ライブラリQCトラブルシューティング

### 9.1 低収量問題

**症状: Total DNA <100 ng**

```yaml
原因分析:

  原因1: 元のDNA/RNA量不足
    確認:
      - cfDNA/cfRNA抽出のQubit値
      - サンプル品質 (溶血、凝固)

    対処:
      - サンプル量増加 (5 mL → 10 mL血漿)
      - 複数抽出バッチのプール
      - 濃縮 (SpeedVac)

  原因2: AMPure XP精製ロス大
    確認:
      - 各精製ステップ前後のQubit測定
      - 回収率計算

    対処:
      - AMPure比率調整 (0.4× → 0.6×)
      - 溶出時間延長 (5分 → 10分)
      - 溶出温度最適化 (50°C)

  原因3: Adapter ligation失敗
    確認:
      - qPCR with adapter primer
      - Ligation効率 <50%

    対処:
      - Ligase fresh試薬使用
      - Ligation時間延長 (10分 → 15分)
      - Adapter量増加 (1.5倍)

  原因4: DNAの分解
    確認:
      - TapeStationでsmear増加
      - 保存期間長い

    対処:
      - 新鮮サンプル使用
      - 凍結融解回数削減
```

### 9.2 N50低下問題

**症状: N50 <5 kb**

```yaml
原因分析:

  原因1: サンプルDNA断片化
    確認:
      - DNA抽出直後のN50測定
      - 元から短い

    対処:
      - サンプル取り扱い改善
      - 抽出プロトコル最適化
      - 短鎖でも解析続行可能

  原因2: AMPure XP不適切
    確認:
      - AMPure前後のN50比較
      - 大幅低下

    対処:
      - AMPure比率変更 (0.4× → 0.6×)
      - 長鎖DNA優先回収

  原因3: ピペッティングによるシェアリング
    確認:
      - 操作後のN50低下

    対処:
      - ピペッティング回数削減
      - Wide bore tipsデータ:
      - ボルテックス厳禁
```

### 9.3 アダプターダイマー過剰

**症状: <200 bp画分 >10%**

```yaml
原因と対処:

  原因1: Adapter過剰
    対処:
      - Adapter量50%に削減
      - cDNA量増加

  原因2: Ligation時間過剰
    対処:
      - 10分 → 5分に短縮

  原因3: AMPure精製不足
    対処:
      - 0.6× → 0.5× (より厳格)
      - SFB洗浄2回実施
      - 追加AMPure精製ラウンド

  緊急対処 (ライブラリ調製後):
    1. BluePippin size selection:
       - 設定: >3 kb回収
       - 回収率: 70-80%
       - ダイマー完全除去

    2. AMPure追加精製:
       - 0.4× ratio
       - 2ラウンド実施
```

----

## 10. ALCOA+準拠記録

### 10.1 ライブラリQC記録統合シート

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MinION ライブラリ品質評価記録
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

実施日: YYYY-MM-DD
実施者: [フルネーム] _____________ 署名: _______

【1. ライブラリ基本情報】
  ライブラリID: ___________________________________
  サンプルID: _____________________________________
  ライブラリタイプ: □ DNA  □ cDNA (RNA由来)
  調製日: YYYY-MM-DD
  保存条件: □ 4°C  □ -20°C  □ -80°C
  凍結融解回数: _______ 回

【2. Qubit定量】
  測定日時: YYYY-MM-DD HH:MM
  Assay: □ dsDNA HS  □ dsDNA BR
  測定者: _____________

  測定値 (3連測定):
    測定1: _______ ng/μL
    測定2: _______ ng/μL
    測定3: _______ ng/μL
    平均: _______ ng/μL
    CV: _______ %

  総量計算:
    体積: _______ μL
    総量: _______ ng

  判定:
    □ 優良 (30-50 ng/μL, >360 ng)
    □ 良好 (20-30 ng/μL, 240-360 ng)
    □ 許容 (10-20 ng/μL, 120-240 ng)
    □ 不合格 (<10 ng/μL, <120 ng)

【3. TapeStation/Bioanalyzer】
  測定日時: YYYY-MM-DD HH:MM
  システム: □ TapeStation D5000  □ Bioanalyzer
  測定者: _____________

  結果:
    Read N50: _______ bp
    平均サイズ: _______ bp
    DIN (DNA Integrity Number): _______

    サイズ分布:
      <500 bp: _______ %
      500 bp - 3 kb: _______ %
      3 kb - 10 kb: _______ %
      10 kb - 20 kb: _______ %
      >20 kb: _______ %

    アダプターダイマー (100-300 bp):
      割合: _______ %
      判定: □ <2% (優良)  □ 2-5% (良好)
            □ 5-10% (許容)  □ >10% (不合格)

  電気泳動イメージ添付: □

  判定:
    □ 優良 (N50 >12 kb, ダイマー <2%)
    □ 良好 (N50 8-12 kb, ダイマー 2-5%)
    □ 許容 (N50 5-8 kb, ダイマー 5-10%)
    □ 不合格 (N50 <5 kb, ダイマー >10%)

【4. NanoDrop純度】
  測定日時: YYYY-MM-DD HH:MM
  測定者: _____________

  結果:
    濃度: _______ ng/μL
    A260/A280: _______
      判定: □ 1.8-2.0 (高純度)
            □ 1.6-1.8 (やや低い)
            □ <1.6 (タンパク質汚染)

    A260/A230: _______
      判定: □ 2.0-2.2 (高純度)
            □ 1.5-2.0 (やや低い)
            □ <1.5 (塩類汚染)

  吸光スペクトル評価:
    □ 正常 (260 nmピーク、smoothcurve)
    □ 異常 (異常ピークあり、記述: ___________)

【5. Molarity計算】
  Qubit濃度: _______ ng/μL
  TapeStation平均サイズ: _______ bp

  計算:
    Molarity = (濃度 × 10^6) / (サイズ × 650)
    Molarity = _______ fmol/μL

  ローディング量計算:
    目標fmol: _______ fmol (通常75 fmol)
    必要体積: _______ μL

    判定: □ 適切 (8-25 μL)
          □ 要希釈 (>25 μL)
          □ 不足 (<8 μL)

【6. Duplex予測評価】
  DNA完全性評価: □ >90%  □ 80-90%  □ <80%
  N50値: _______ kb
  予測Duplex rate: _______ %

  判定:
    □ 優良 (50-60%期待)
    □ 良好 (40-50%期待)
    □ 許容 (30-40%期待)
    □ 低い (<30%期待)

【7. 総合判定】
  □ フローセルローディング可 (優良品質)
  □ ローディング可 (標準品質)
  □ 条件付き可 (理由: ______________________)
  □ 不可、要改善 (理由: ____________________)

  改善必要項目:
    □ アダプターダイマー除去 (AMPure追加)
    □ 濃度不足 (複数バッチプール)
    □ N50低下 (サイズ選択)
    □ 純度低下 (追加精製)
    □ ライブラリ再調製

【8. 次のステップ】
  □ フローセルプライミングへ進む
  □ ライブラリ希釈実施
  □ AMPure追加精製
  □ サイズ選択 (BluePippin)
  □ ライブラリ再調製
  □ 保存 (温度: ____ °C、期間: ____ 日)

【9. 特記事項】
  ________________________________________________
  ________________________________________________
  ________________________________________________

【10. 確認・承認】
  測定者: [フルネーム] _____________ 署名: _______
  日時: YYYY-MM-DD HH:MM

  確認者: [フルネーム] _____________ 署名: _______
  日時: YYYY-MM-DD HH:MM

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

----

## まとめ

本章では、MinIONライブラリの品質評価とDuplex対応について詳述しました。

**重要ポイント:**

1. **ライブラリ品質の重要性**
   - N50 >10 kb、アダプターダイマー <5%
   - Duplex率50%達成の鍵
   - PMDA病原体検出精度に直結

2. **Duplex Sequencingの必須性**
   - Q30+データで偽陽性劇的削減
   - SNP level病原体同定可能
   - 系統解析、薬剤耐性検出

3. **多層QCシステム**
   - Qubit (定量)
   - TapeStation (サイズ分布)
   - NanoDrop (純度)
   - 総合評価で判定

4. **Duplex最適化戦略**
   - DNA完全性維持
   - 長鎖DNA回収
   - 高効率Adapter ligation
   - 48時間シーケンスラン

5. **品質管理の標準化**
   - バッチ間変動 <15%
   - Positive/Negative control
   - 内部標準使用

6. **ALCOA+準拠記録**
   - 全QCデータ記録
   - トレーサビリティ確保
   - 長期データ保存

**次章予告:**

第9章では、重要QC基準と合否判定について詳述します。各ステップの合格基準、逸脱時の対処、QC failの根本原因分析について解説します。

----

**文書情報**
- 作成日: 2025-03-08
- バージョン: 1.0
- 作成者: MinIONメタゲノム解析プロトコル開発チーム
- 承認者: [承認者名]
- 次回改訂予定: 2025-09-08
