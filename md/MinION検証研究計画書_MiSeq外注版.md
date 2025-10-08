# MinION病原体検査法検証研究計画書
## MiSeq外注比較による最適プロトコル確立とPMDA承認体制構築（改訂版）

**研究題目**: 異種移植用ドナーブタにおけるMinION単独病原体メタゲノム解析法の分析的バリデーション
**研究期間**: 12ヶ月（Month 1-12）
**研究デザイン**: 前向き比較検証研究（MinION内製 vs MiSeq外注並行解析）
**最終目標**: MinION単独によるPMDA承認可能な病原体検査体制の完全内製化

**作成日**: 2025年10月8日
**バージョン**: 2.0（MiSeq外注対応版）

----

## エグゼクティブサマリー

### 研究概要（改訂版）

```yaml
【研究の位置づけ】

現状:
  - MinION: 新規技術、臨床検査実績少ない
  - MiSeq: ゴールドスタンダード、PMDA承認実績豊富
  - 戦略: MiSeq解析は外注、MinIONのみ内製化

研究戦略:
  Phase_1_並行検証（Month 1-7）:
    - MinION解析: 自社内製
    - MiSeq解析: 外部受託機関に委託
    - MiSeqを「Reference Standard」として使用
    - MinIONプロトコル・解析パイプライン最適化
    - 一致率・精度評価

  Phase_2_独立検証（Month 8-10）:
    - 最適化MinIONプロトコルの独立検証
    - MiSeq外注結果と盲検比較
    - 統計的同等性証明

  Phase_3_完全内製化（Month 11-12）:
    - MinION単独運用開始
    - MiSeq外注は品質管理のみ（四半期1回）
    - PMDA提出資料作成

Exit Criteria（MiSeq外注終了の判断基準）:
  1. ✓ 種レベル同定一致率 >95%
  2. ✓ 定量相関係数 R² >0.90
  3. ✓ 陽性一致率（PPA） >95%
  4. ✓ 陰性一致率（NPA） >98%
  5. ✓ MinION独立検証で全基準達成

最終成果:
  【MinION完全内製化達成】
    - MiSeq装置不要（初期投資削減: -15M円）
    - 初期投資: 8.4M円（MinION本体のみ）
    - 年間運営: 18.6M円（外注費最小化後）
    - PMDAバリデーションデータ完備
    - SOP・解析パイプライン確立
    - Phase 3以降はMiSeq外注ほぼゼロ

投資削減効果:
  従来計画（MiSeq購入）: 初期49.4M円
  改訂計画（MiSeq外注）: 初期17.4M円
  削減額: -32M円（65%削減）
```

----

## 1. 研究目的と仮説

### 1.1 主要目的（Primary Objective）

```yaml
Primary_Objective:

  目的:
    "MinION内製解析が、MiSeq外注解析（Reference Standard）と
     同等の精度でPMDA指定91病原体を検出・定量できることを証明する"

  主要評価項目:
    1. 種レベル同定一致率（Species-level concordance）
    2. 定量的一致性（Quantitative agreement）

  成功基準:
    - 一致率 >95%（95% CI下限 >90%）
    - 定量相関 R² >0.90
    - Bland-Altman plot: Bias <±0.3 log₁₀

  意義:
    "外注MiSeqとの同等性証明により、
     MinION内製のみでPMDA承認取得可能となる"
```

### 1.2 副次目的（Secondary Objectives）

```yaml
Secondary_Objectives:

  1_感度・特異度評価:
    目的: "MinION内製解析の診断性能評価"
    評価項目:
      - Positive Percent Agreement (PPA)
      - Negative Percent Agreement (NPA)
      - LOD（Limit of Detection）

  2_最適プロトコル確立:
    目的: "MinION性能を最大化する内製プロトコル決定"
    比較条件:
      - Basecalling mode（Simplex vs Duplex）
      - Sequencing depth（20Gb vs 30Gb vs 50Gb）
      - Library prep kit（Rapid vs Ligation）

  3_内製化のコスト効率:
    目的: "完全内製化の経済性評価"
    評価項目:
      - Per sample cost（MinION内製 vs MiSeq外注）
      - Phase 3以降の年間総コスト
      - 外注費削減額

  4_ターンアラウンドタイム優位性:
    目的: "MinION内製の迅速性実証"
    評価項目:
      - MinION内製 TAT（目標 <48時間）
      - MiSeq外注 TAT（実測値記録）
      - 差分の臨床的価値

  5_PERV詳細解析（MinION独自価値）:
    目的: "内製化による高度解析能力獲得"
    評価項目:
      - 全ゲノム取得率
      - 組換え部位同定精度
      - MiSeq外注では得られない情報の価値
```

----

## 2. 研究デザイン（改訂版）

### 2.1 全体スキーム

```yaml
Study_Design:

  Type: "Prospective Comparative Validation Study"

  Comparison: "In-house MinION vs Outsourced MiSeq"

  Key_Change: "MiSeq解析を外部受託機関に委託"

  Blinding:
    Phase_1（最適化期）: "非盲検"
    Phase_2（独立検証期）: "解析者盲検化"

  Duration: "12ヶ月"

  Sample_size: "60サンプル"

  #========== 研究フロー（改訂版） ==========
  Study_Flow:

    Phase_1_Protocol_Optimization（Month 1-7）:
      目的: "MinION最適条件決定"
      サンプル数: "30サンプル"

      Workflow:
        1: "同一サンプルを2分割"
           ├─ MinION用: 自社で解析（3条件）
           └─ MiSeq用: 外部受託機関に送付

        2: "MinION 3条件を自社で比較"
           ├─ Condition A: Simplex, 20Gb, Rapid kit
           ├─ Condition B: Duplex, 30Gb, Ligation kit
           └─ Condition C: Duplex, 50Gb, Ligation kit

        3: "MiSeq外注結果を正解として一致率算出"

        4: "最適条件決定（Month 7）"

      外注先との連携:
        - サンプル送付（週1-2回バッチ）
        - 標準化プロトコル指定
        - 生データ（FASTQ）受領
        - 結果レポート受領

      分析:
        - MinION 3条件間比較
        - 各条件とMiSeq外注結果の一致率
        - コスト-性能バランス評価

      成果物:
        - 最適化MinIONプロトコル（SOP）
        - 自社解析パイプライン確定版

    Phase_2_Independent_Validation（Month 8-10）:
      目的: "最適MinION内製解析の独立検証"
      サンプル数: "30サンプル（新規）"

      Workflow:
        1: "Phase 1最適条件でMinION自社解析"
        2: "並行してMiSeq外注"
        3: "解析者盲検化"
        4: "完了後に結果突合"

      判定基準:
        Success: "全Exit Criteria達成 → Phase 3移行"
        Failure: "Phase 1に戻り再最適化"

    Phase_3_Full_In_House_Operation（Month 11-12）:
      目的: "MinION完全内製運用開始"
      サンプル数: "臨床サンプル（実運用）"

      Workflow:
        1: "MinION内製のみで解析"
        2: "MiSeq外注は以下のみ:"
           - QC目的（四半期1サンプル）
           - MinION陽性の確認（疑義例のみ、年1-2回）

      成果:
        - バリデーション報告書
        - 内製SOP最終版
        - PMDA相談資料
        - 外注依存からの完全脱却
```

### 2.2 MiSeq外注先の選定基準

```yaml
External_Lab_Selection:

  #===== 必須要件 =====
  Mandatory_Requirements:

    1_技術要件:
      - "Illumina MiSeq保有"
      - "メタゲノム解析実績（病原体検出）"
      - "FASTQ生データ提供可能"
      - "標準化プロトコル対応"

    2_品質要件:
      - "ISO 15189認定 or 同等の品質管理"
      - "データ品質保証（Q30 >80%）"
      - "トレーサビリティ確保"

    3_セキュリティ要件:
      - "秘密保持契約（NDA）締結"
      - "データ暗号化転送"
      - "保管期間・破棄規定明確"

    4_納期要件:
      - "TAT 2-3週間（サンプル受領から結果報告）"
      - "緊急対応可能（オプション）"

  #===== 推奨外注先 =====
  Recommended_Vendors:

    Vendor_A: "Macrogen Japan"
      所在地: "東京"
      実績: "国内最大手、メタゲノム実績豊富"
      TAT: "2-3週間"
      価格: "約150,000円/サンプル"
      利点: "信頼性高、PMDA提出実績あり"

    Vendor_B: "タカラバイオ"
      所在地: "滋賀・東京"
      実績: "国内メーカー、カスタム対応"
      TAT: "2-3週間"
      価格: "約180,000円/サンプル"
      利点: "日本企業、きめ細かい対応"

    Vendor_C: "GENEWIZ Japan"
      所在地: "東京"
      実績: "グローバル企業、標準化"
      TAT: "2-3週間"
      価格: "約120,000円/サンプル"
      利点: "コスト競争力"

  #===== 選定プロセス =====
  Selection_Process:

    Step_1_RFP発行:
      - "要求仕様書作成（技術・品質・納期）"
      - "3社以上に見積依頼"

    Step_2_評価:
      配点:
        技術力: "30点"
        品質管理: "30点"
        価格: "20点"
        納期: "10点"
        実績: "10点"

    Step_3_契約:
      - "最高得点の1社と契約"
      - "バックアップとして次点と基本合意"

  #===== 発注仕様 =====
  Order_Specification:

    プロトコル指定:
      Platform: "MiSeq"
      Reagent: "MiSeq Reagent Kit v3 (600 cycles)"
      Read_config: "2 × 300bp Paired-end"
      Depth: "2M reads/sample minimum"
      PhiX: "1% spike-in"

    納品物:
      1: "FASTQ files（Raw data）"
      2: "QCレポート（FastQC, MultiQC）"
      3: "解析結果"
         - Host removal済みリード
         - 病原体検出結果（Kraken2/BLAST）
         - 定量結果（copies/mL）
      4: "最終レポート（PDF）"

    データ転送:
      方法: "専用FTPサーバー（暗号化）"
      容量: "5-10Gb/サンプル"
```

----

## 3. 実験プロトコル（改訂版）

### 3.1 サンプル前処理（自社）

```yaml
Sample_Processing:

  目的: "MinION内製用とMiSeq外注用に分割"

  #===== Step 1: 血漿分離 =====
  Plasma_Separation:

    採血:
      容量: "20mL"
      抗凝固剤: "EDTA-2K"
      保存: "4°C、24時間以内に処理"

    遠心分離:
      条件: "1,500×g、10分、4°C"
      上清回収: "15mL血漿"

    分注（重要な変更）:
      MinION_in_house用: "10mL"
      MiSeq_outsource用: "4mL"
      予備: "1mL"

  #===== Step 2: cfDNA/RNA抽出 =====
  Nucleic_Acid_Extraction:

    Kit: "Zymo Quick-cfDNA/RNA Serum & Plasma Kit"

    Input:
      MinION用: "10mL血漿 → 自社で抽出"
      MiSeq用: "4mL血漿 → 自社で抽出"

    Protocol: "製造者推奨プロトコル"

    Expected_yield:
      MinION用: "30-100ng"
      MiSeq用: "10-40ng"

    QC:
      Qubit測定: "両サンプル"
      TapeStation: "両サンプル"

      合格基準:
        - "Total NA >20ng（MinION）"
        - "Total NA >10ng（MiSeq）"
        - "Mean size >500bp"

  #===== Step 3: ホスト除去 =====
  Host_Depletion:

    Kit: "NEBNext Microbiome DNA Enrichment Kit"

    Input:
      MinION用: "抽出cfDNA/RNA全量"
      MiSeq用: "抽出cfDNA/RNA全量"

    Expected_depletion: "90-99%"

    QC:
      qPCR: "Sus scrofa β-actin定量"
      目標: ">95% depletion"

    Output:
      MinION用: "病原体濃縮NA（10-20ng） → 即座にライブラリ調製"
      MiSeq用: "病原体濃縮NA（5-10ng） → 凍結保存 → バッチ送付"
```

### 3.2 MinION内製プロトコル（自社）

```yaml
MinION_In_House_Protocol:

  Location: "自社ラボ"

  #===== 共通設定 =====
  Common_Settings:

    Device: "MinION Mk1D"
    Flowcell: "R10.4.1"
    Chemistry: "Kit 14"

  #========== Condition A: 高速・低コスト ==========
  Condition_A_Rapid:

    Library_Prep:
      Kit: "Rapid Sequencing Kit (SQK-RAD114)"
      Input: "200ng"
      Time: "10分"

    Sequencing:
      Duration: "6時間"
      Expected_output: "20Gb"

    Basecalling:
      Location: "AWS EC2 G5 GPU"
      Mode: "Simplex"
      Model: "dna_r10.4.1_e8.2_400bps_fast"
      Expected_Q: "Q15-18"

    Analysis:
      Location: "AWS（自社パイプライン）"
      Pipeline: "Nextflow"
      Steps:
        - QC
        - Host removal（Sus scrofa）
        - Pathogen detection（Kraken2/BLAST）
        - Quantification
        - Report generation

    Total_TAT: "24-36時間（前処理 → レポート）"

    Cost_per_sample: "約11,000円（試薬 + AWS）"

  #========== Condition B: バランス型（推奨候補） ==========
  Condition_B_Standard:

    Library_Prep:
      Kit: "Ligation Sequencing Kit (SQK-LSK114)"
      Input: "200ng"
      Time: "2時間"

    Sequencing:
      Duration: "24時間"
      Expected_output: "30Gb"

    Basecalling:
      Location: "AWS EC2 G5 GPU"
      Mode: "Duplex（高精度）"
      Model: "dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
      Expected_Q: "Q28-30"

    Total_TAT: "48-60時間"

    Cost_per_sample: "約28,000円（試薬 + AWS）"

  #========== Condition C: 最高精度 ==========
  Condition_C_High_Accuracy:

    Sequencing:
      Duration: "48時間"
      Expected_output: "50Gb"

    Total_TAT: "72-84時間"

    Cost_per_sample: "約43,000円（試薬 + AWS）"
```

### 3.3 MiSeq外注プロトコル

```yaml
MiSeq_Outsourcing_Protocol:

  Location: "外部受託機関"

  #===== サンプル送付準備（自社） =====
  Sample_Shipment_Preparation:

    Step_1_凍結保存:
      - "Host depletion後のNA"
      - "-80°C凍結"
      - "ドライアイス送付まで保管"

    Step_2_バッチング:
      頻度: "週1回 or 5サンプル貯まったら"
      理由: "送料削減、効率化"

    Step_3_梱包:
      容器: "1.5mL凍結チューブ"
      梱包: "ドライアイス、断熱ボックス"
      ラベル: "サンプルID、濃度、日付"

    Step_4_書類:
      - "検査依頼書（サンプルリスト）"
      - "安全データシート（該当する場合）"
      - "送付状"

    Step_5_輸送:
      方法: "宅配便（ヤマト運輸 クール便）"
      時間: "翌日配達"
      追跡: "配送番号記録"

  #===== 外注先での解析 =====
  External_Lab_Analysis:

    受領確認:
      - "外注先からサンプル到着確認メール"
      - "数量・品質確認"

    解析実施:
      Protocol: "事前合意の標準プロトコル"

      Library_Prep:
        Kit: "Nextera XT DNA Library Prep Kit"
        Input: "1ng"

      Sequencing:
        Platform: "MiSeq"
        Reagent: "v3 600 cycles"
        Read: "2 × 300bp"
        Depth: ">2M reads/sample"

      Data_Analysis:
        Pipeline: "外注先標準パイプライン"
        Steps:
          - QC（FastQC）
          - Host removal（BWA）
          - Pathogen detection（Kraken2/BLAST）
          - Quantification

    TAT:
      サンプル受領: "Day 0"
      解析完了: "Day 14-21"
      結果報告: "Day 21"

  #===== データ受領（自社） =====
  Data_Reception:

    納品物:
      1_Raw_data:
        - "FASTQ files（暗号化FTP）"
        - "サイズ: 5-10Gb/サンプル"

      2_QC_report:
        - "FastQC/MultiQC結果"
        - "Run metrics"

      3_Analysis_results:
        - "病原体検出リスト（TSV）"
        - "定量結果（copies/mL）"

      4_Final_report:
        - "PDFレポート"

    受領確認:
      - "データ完全性チェック（MD5 checksum）"
      - "ファイル数・サイズ確認"
      - "レポート内容確認"

    保存:
      - "AWS S3にバックアップ"
      - "ローカルNASにも保存（二重化）"
```

### 3.4 並行解析の実施手順（改訂版）

```yaml
Parallel_Analysis_SOP:

  #===== サンプル処理タイムライン =====
  Sample_Processing_Timeline:

    Day_0（自社ラボ）:
      - "採血"
      - "血漿分離"
      - "2分割（MinION用 + MiSeq用）"
      - "cfDNA/RNA抽出"
      - "Host depletion"
      - "QC確認"

    Day_1（自社ラボ）:
      MinION_pathway:
        - "ライブラリ調製開始"
        - "シーケンシング開始"

      MiSeq_pathway:
        - "-80°C凍結保存"
        - "バッチ待機"

    Day_1-3（MinION自社解析）:
      - "シーケンシング継続"
      - "AWS basecalling"
      - "解析パイプライン実行"
      - "Day 3: MinION結果完成"

    Day_7（MiSeq外注準備）:
      - "5サンプル貯まった時点"
      - "バッチ梱包"
      - "外注先へ発送"

    Day_8-28（MiSeq外注期間）:
      - "Day 8: 外注先到着"
      - "Day 8-28: 解析実施"
      - "Day 28: 結果受領"

    Day_29（データ統合）:
      - "MinION vs MiSeq結果比較"
      - "一致率算出"
      - "データベース登録"

  #===== データ管理 =====
  Data_Management:

    LIMS（自社）:
      Tool: "Excel or 簡易データベース"

      記録項目:
        - "サンプルID"
        - "採取日時"
        - "分注日"
        - "MinION解析日"
        - "MiSeq発送日"
        - "MiSeq結果受領日"
        - "外注先名"
        - "QC結果"

    Raw_Data:
      MinION: "AWS S3 + ローカルNAS"
      MiSeq: "外注先FTP → AWS S3ダウンロード"
      保存期間: "5年間"

    Results:
      保存先: "AWS RDS PostgreSQL + Excel"
      バックアップ: "週次"
```

----

## 4. データ解析計画（改訂版）

### 4.1 主要評価項目の解析

```yaml
Primary_Endpoint_Analysis:

  #===== Endpoint 1: 種レベル同定一致率 =====
  Species_Concordance:

    比較対象:
      Reference: "MiSeq外注結果"
      Test: "MinION内製結果"

    定義:
      "MinION内製とMiSeq外注で同一の病原体種を検出した割合"

    計算方法: |
      Concordance (%) = (一致検出数) / (総比較数) × 100

      一致の定義:
        - 両方とも同じ種を検出（True Positive）
        - 両方とも陰性（True Negative）

      不一致の分類:
        - MinION陽性、MiSeq陰性（False Positive）
        - MinION陰性、MiSeq陽性（False Negative）

    統計解析:
      - "一致率の95%信頼区間（Wilson score）"
      - "McNemar検定"
      - "Cohen's Kappa係数"

    成功基準:
      - "一致率 >95%（95% CI下限 >90%）"
      - "Kappa >0.90"

  #===== Endpoint 2: 定量的一致性 =====
  Quantitative_Agreement:

    対象: "両プラットフォームで陽性となった病原体"

    解析手法:

      1_相関分析:
        X軸: "MiSeq外注定量値（log₁₀ copies/mL）"
        Y軸: "MinION内製定量値（log₁₀ copies/mL）"
        Method: "Pearson/Spearman相関"
        Plot: "Scatter plot with regression line"
        成功基準: "R² >0.90"

      2_Bland_Altman_plot:
        X軸: "(MinION + MiSeq) / 2"
        Y軸: "MinION - MiSeq"
        計算: "Mean bias, 95% LoA"

        成功基準:
          - "Mean bias: -0.3 ~ +0.3 log₁₀"
          - "95% LoA: ±1.0 log₁₀以内"

      3_Passing_Bablok回帰:
        成功基準:
          - "傾き: 0.9-1.1（95% CI）"
          - "切片: -0.3 ~ +0.3（95% CI）"
```

### 4.2 副次評価項目の解析

```yaml
Secondary_Endpoint_Analysis:

  #===== PPA/NPA（診断性能） =====
  Diagnostic_Performance:

    Reference_standard: "MiSeq外注結果"

    定義:

      PPA（Positive Percent Agreement）:
        計算: "PPA = TP / (TP + FN) × 100"
        意味: "MiSeq陽性をMinIONが検出できる割合"

      NPA（Negative Percent Agreement）:
        計算: "NPA = TN / (TN + FP) × 100"
        意味: "MiSeq陰性をMinIONも陰性とする割合"

    成功基準:
      - "PPA >95%"
      - "NPA >98%"

  #===== コスト比較 =====
  Cost_Analysis:

    MinION_内製:
      Phase_1（3条件平均）:
        試薬: "27,000円/サンプル"
        AWS: "3,000円/サンプル"
        人件費: "5,000円/サンプル"
        Total: "35,000円/サンプル"

      Phase_2以降（最適条件）:
        試薬: "25,000円/サンプル"
        AWS: "3,000円/サンプル"
        人件費: "3,000円/サンプル（習熟）"
        Total: "31,000円/サンプル"

    MiSeq_外注:
      Phase_1-2:
        外注費: "150,000円/サンプル"
        前処理: "8,000円/サンプル（自社）"
        送料: "2,000円/サンプル"
        Total: "160,000円/サンプル"

    削減効果:
      Per_sample: "160,000 - 31,000 = 129,000円/サンプル"
      削減率: "81%"

      Phase_3以降（年24サンプル）:
        MinION内製: "31,000 × 24 = 744,000円/年"
        MiSeq外注（全量）: "160,000 × 24 = 3,840,000円/年"
        削減額: "3,096,000円/年"

  #===== TAT比較 =====
  Turnaround_Time_Comparison:

    MinION_内製:
      Condition_A: "24-36時間"
      Condition_B: "48-60時間"
      Condition_C: "72-84時間"

    MiSeq_外注:
      実測値記録:
        前処理（自社）: "1日"
        発送: "1日"
        外注先解析: "14-21日"
        結果受領: "1日"
        Total: "17-24日"

    優位性:
      MinION勝利: "15-23日短縮"
      臨床的価値: "移植判断の迅速化"
```

----

## 5. MinIONプロトコル最適化戦略（改訂版）

### 5.1 Phase 1最適化プロセス

```yaml
Optimization_Process:

  #===== 評価軸（Multi-criteria decision analysis） =====
  Evaluation_Criteria:

    1_精度（Weight: 50%）:
      指標:
        - "MiSeq外注との種同定一致率"
        - "定量相関R²"
        - "PPA/NPA"
      目標: "一致率 >95%"
      重要性: "最重要（PMDA承認要件）"

    2_コスト（Weight: 30%）:
      指標: "Per sample cost（試薬+AWS+人件費）"
      目標: "MiSeq外注（160,000円）より大幅に安価"
      重要性: "内製化の主目的"

    3_TAT（Weight: 15%）:
      指標: "Total time（受領→レポート）"
      目標: "<72時間"
      重要性: "緊急対応能力"

    4_操作性（Weight: 5%）:
      指標:
        - "Hands-on time"
        - "失敗率"
      目標: "再現性 >95%"

  #===== スコアリング =====
  Scoring_System:

    Condition_A（Rapid, Simplex, 20Gb）:
      精度: "85点（一致率92%想定）"
      コスト: "100点（11,000円 vs 外注160,000円）"
      TAT: "100点（24-36時間）"
      操作性: "100点（10分調製）"

      Total: "85×0.5 + 100×0.3 + 100×0.15 + 100×0.05 = 92.5点"

      判定: "コスト・TAT優秀、精度やや不足"

    Condition_B（Ligation, Duplex, 30Gb）:
      精度: "98点（一致率97%想定）"
      コスト: "95点（28,000円 vs 外注160,000円）"
      TAT: "90点（48-60時間）"
      操作性: "90点（2時間調製）"

      Total: "98×0.5 + 95×0.3 + 90×0.15 + 90×0.05 = 95.5点"

      判定: "✓ バランス最高、推奨候補"

    Condition_C（Ligation, Duplex, 50Gb）:
      精度: "100点（一致率98%想定）"
      コスト: "85点（43,000円）"
      TAT: "80点（72-84時間）"
      操作性: "90点（2時間調製）"

      Total: "100×0.5 + 85×0.3 + 80×0.15 + 90×0.05 = 91.5点"

      判定: "精度最高、コスト・TATやや劣る"

  #===== 予想される最適条件 =====
  Expected_Optimal:

    第1候補: "Condition B（バランス型）"
    理由:
      - "精度十分（97%、MiSeq外注と同等）"
      - "コスト: 外注の82%削減"
      - "TAT: 外注の1/7"
      - "総合スコア最高"

    最終判定: "Phase 1実測データで確定"
```

----

## 6. 品質管理・品質保証（改訂版）

### 6.1 QC/QAシステム

```yaml
QC_QA_System:

  #===== MinION内製QC =====
  MinION_In_House_QC:

    Controls:

      Positive_Control:
        頻度: "各ラン必須"
        内容: "既知病原体DNA（PERV, E.coli）"
        濃度: "500 copies/mL"
        判定: "検出必須、定量±30%以内"

      Negative_Control:
        頻度: "各ラン必須"
        内容: "DNase/RNase処理水"
        判定: "病原体検出ゼロ"

      Internal_Control（Spike-in）:
        頻度: "全サンプル"
        内容: "Lambda phage DNA"
        濃度: "1,000 copies/mL"
        判定: "回収率 80-120%"

    QC_Metrics:
      Pre-sequencing:
        - "Flowcell active pore数 >1,200"
        - "Library濃度 >10ng/μL"

      Post-sequencing:
        - "Total bases >閾値"
        - "Mean Q >閾値"
        - "Spike-in回収率合格"

  #===== MiSeq外注QC =====
  MiSeq_Outsourcing_QC:

    発注時チェック:
      - "サンプル濃度 >5ng"
      - "純度確認（260/280比）"
      - "ラベル確認"

    外注先QCレポート確認:
      必須項目:
        - "Total reads >2M"
        - "Q30 >80%"
        - "PhiX error rate <1%"

      不合格時:
        - "外注先に再解析依頼"
        - "追加料金なしで再実施（契約条項）"

    データ受領時チェック:
      - "ファイル完全性（checksum）"
      - "FASTQ形式確認"
      - "サンプルID一致確認"

  #===== 比較QC（MinION vs MiSeq外注） =====
  Comparative_QC:

    Phase_1-2期間:

      月次レビュー:
        確認項目:
          - "一致率トレンド"
          - "不一致サンプルの分析"
          - "外注先品質安定性"

      四半期監査:
        - "外注先訪問（可能であれば）"
        - "プロトコル遵守確認"
        - "QC記録閲覧"

    Phase_3以降:

      外注頻度: "四半期1サンプル"
      目的: "MinION内製の継続的品質確認"

      判定基準:
        - "一致率維持 >95%"
        - "不一致時は原因調査 → CAPA"
```

### 6.2 外注先管理

```yaml
Vendor_Management:

  #===== 契約管理 =====
  Contract_Management:

    契約期間: "1年間（自動更新）"

    契約内容:
      納期: "サンプル受領後21日以内"
      品質保証: "Q30 >80%保証、不合格時再解析無料"
      秘密保持: "NDA、データ破棄規定"
      価格: "150,000円/サンプル（Phase 1-2: 60サンプル）"

    SLA（Service Level Agreement）:
      納期遵守率: ">90%"
      品質合格率: ">95%"
      未達時: "ペナルティ or 値引き"

  #===== 定期レビュー =====
  Periodic_Review:

    頻度: "四半期ごと"

    評価項目:
      - "納期遵守率"
      - "品質合格率"
      - "MinION一致率"
      - "コミュニケーション品質"

    Action:
      良好: "契約継続"
      問題あり: "改善要求"
      重大問題: "契約解除、別業者切替"

  #===== バックアップ体制 =====
  Backup_Vendor:

    準備:
      - "次点業者と基本合意"
      - "緊急時3営業日以内に切替可能"

    発動条件:
      - "主契約先の重大品質問題"
      - "納期遅延常習（3回以上）"
      - "倒産・廃業"
```

----

## 7. タイムライン・マイルストーン（改訂版）

### 7.1 詳細スケジュール

```yaml
Detailed_Timeline:

  #========== Phase 1: Protocol Optimization ==========
  Phase_1_Month_1_to_7:

    Month_1:
      Week_1_2:
        □ "研究計画書最終化"
        □ "外注先選定（RFP発行、評価、契約）"
        □ "MinION調達完了確認"
        □ "試薬発注"
        □ "AWS環境構築"

      Week_3_4:
        □ "外注先とプロトコル合意"
        □ "初回サンプル5個でパイロット試験"
        □ "MinION解析パイプライン動作確認"
        □ "外注先へ初回サンプル送付テスト"

    Month_2:
      □ "Phase 1サンプル収集開始（n=30）"
      □ "週2-3サンプルペース"
      □ "MinION 3条件解析（自社）"
      □ "MiSeq外注（週1回バッチ送付）"

    Month_3-4:
      □ "継続サンプル解析"
      □ "MinION結果: 1週間以内"
      □ "MiSeq結果: 3-4週間後受領"

    Month_5:
      □ "Phase 1サンプル30個完了"
      □ "全MinION解析完了"
      □ "MiSeq外注結果待ち（最終バッチ）"

    Month_6:
      □ "全MiSeq外注結果受領完了"
      □ "データ統合"

    Month_7:
      □ "Phase 1統計解析"
      □ "最適MinION条件決定"
      □ "SOP草案完成"
      □ "Phase 1報告書作成"
      □ "Go/No-Go判断（Phase 2移行）"

  #========== Phase 2: Independent Validation ==========
  Phase_2_Month_8_to_10:

    Month_8:
      Week_1:
        □ "Phase 2プロトコル最終確認"
        □ "解析者盲検化手順確立"
        □ "新規サンプル収集開始"

      Week_2_4:
        □ "サンプル解析開始（最適MinION条件のみ）"
        □ "MiSeq外注並行"
        □ "週3-4サンプルペース"

    Month_9:
      □ "継続サンプル解析"
      □ "30サンプル完了"
      □ "MinION結果完了"
      □ "MiSeq外注結果待ち"

    Month_10:
      Week_1_2:
        □ "全MiSeq外注結果受領"
        □ "盲検解除"
        □ "MinION vs MiSeq結果突合"

      Week_3:
        □ "統計解析実施"
        □ "Exit criteria達成確認"

      Week_4:
        □ "Phase 2報告書作成"
        □ "Phase 3移行判断"

  #========== Phase 3: Full In-House Operation ==========
  Phase_3_Month_11_to_12:

    Month_11:
      Week_1_2:
        □ "統合バリデーション報告書作成"
        □ "PMDA相談資料準備"

      Week_3_4:
        □ "MinION単独運用SOP最終化"
        □ "MiSeq外注縮小計画策定（四半期1回）"
        □ "トレーニングプログラム実施"

    Month_12:
      □ "MinION完全内製運用開始"
      □ "MiSeq外注なし（この期間）"
      □ "実運用データ収集"
      □ "年次総括報告書作成"
      □ "PMDA提出資料完成"
      □ "次年度計画策定"

  #===== マイルストーン =====
  Key_Milestones:

    M1（Month 1）: "外注先契約完了、研究開始"
    M2（Month 7）: "Phase 1完了、最適条件決定"
    M3（Month 10）: "Phase 2完了、MinION性能証明"
    M4（Month 11）: "バリデーション報告書完成"
    M5（Month 12）: "MinION完全内製運用開始"
    M6（Month 12）: "PMDA提出資料完成"
```

### 7.2 リスク管理（改訂版）

```yaml
Risk_Management:

  #===== 主要リスク =====
  Major_Risks:

    Risk_1_外注先の納期遅延:
      確率: "Medium（30%）"
      影響: "Medium（研究遅延1-2ヶ月）"

      軽減策:
        - "契約でSLA明記（納期21日）"
        - "定期的進捗確認"
        - "バックアップ外注先確保"

      対応計画:
        - "遅延3回 → 外注先切替"
        - "期間延長（+2ヶ月）"

    Risk_2_外注先の品質問題:
      確率: "Low（10%）"
      影響: "High（データ信頼性低下）"

      軽減策:
        - "実績ある大手外注先選定"
        - "パイロット試験で品質確認"
        - "QCレポート厳格チェック"

      対応計画:
        - "不合格サンプル再解析要求（無料）"
        - "重大問題時は外注先変更"

    Risk_3_MinIONがExit criteria未達成:
      確率: "Medium（30%）"
      影響: "High（完全内製化不可）"

      軽減策:
        - "Phase 1で十分最適化"
        - "3条件で広範囲カバー"

      対応計画:
        Option_A: "Phase 1再実施（条件追加）"
        Option_B: "内製+外注ハイブリッド継続"
        Option_C: "MinION補助用途に限定"

    Risk_4_サンプル輸送トラブル:
      確率: "Low（5%）"
      影響: "Low（1サンプル損失）"

      軽減策:
        - "ドライアイス十分量使用"
        - "断熱ボックス使用"
        - "翌日配達指定"
        - "予備サンプル保管"

      対応計画:
        - "予備サンプルで再送"

  #===== 新規リスク（外注特有） =====
  Outsourcing_Specific_Risks:

    Risk_5_秘密保持違反:
      確率: "Very Low（<5%）"
      影響: "Critical（知的財産流出）"

      軽減策:
        - "NDA締結（厳格な条項）"
        - "データ破棄確認（証明書）"
        - "信頼できる大手外注先選定"

    Risk_6_外注費用増加:
      確率: "Low（10%）"
      影響: "Medium（予算超過）"

      軽減策:
        - "固定価格契約（60サンプル一括）"
        - "予備費確保（+10%）"

      対応計画:
        - "Phase 2サンプル数削減（30 → 25）"
```

----

## 8. 予算計画（改訂版）

### 8.1 研究費用詳細

```yaml
Budget_Plan:

  #===== 初期投資（Year 0） =====
  Initial_Investment:

    MinION_System:
      MinION_Mk1D_Starter_Pack: "750,000円（本体+フローセル5個）"
      追加Flowcells_15個: "1,350,000円（90,000円×15）"
      Flowcells合計: "20個"
      Sequencing_kits: "1,000,000円"
      小計: "3,100,000円"

    MiSeq_System:
      装置購入: "0円（外注のため不要）"
      削減額: "-15,000,000円"

    共通機器:
      cfDNA抽出装置: "500,000円"
      Host_depletion_kit: "300,000円"
      QC機器（Qubit, TapeStation）: "2,000,000円"
      遠心機・ピペット等: "500,000円"
      小計: "3,300,000円"

    AWS環境:
      初期構築: "0円（従量課金のみ）"

    初期投資合計: "6,400,000円"
    （従来計画24.45M円 → 削減額: -18.05M円、74%削減）

  #===== Phase別運営費 =====
  Operating_Costs:

    Phase_1（Month 1-7、30サンプル）:

      MinION_内製:
        Flowcells: "30×3条件÷12 = 8個 × 90,000円 = 720,000円"
        Library_kits:
          Rapid: "10 × 8,000円 = 80,000円"
          Ligation: "20 × 25,000円 = 500,000円"
        AWS: "3,000円 × 30 = 90,000円"
        小計: "1,390,000円"

      MiSeq_外注:
        外注費: "30サンプル × 150,000円 = 4,500,000円"
        送料: "30 × 2,000円 = 60,000円"
        小計: "4,560,000円"

      共通試薬（自社前処理）:
        cfDNA抽出: "30 × 8,000円 = 240,000円"
        Host_depletion: "30 × 8,000円 = 240,000円"
        QC試薬: "30 × 3,000円 = 90,000円"
        小計: "570,000円"

      Phase_1合計: "6,520,000円"

    Phase_2（Month 8-10、30サンプル）:

      MinION_内製（最適条件のみ）:
        Flowcells: "30÷12 = 3個 × 90,000円 = 270,000円"
        Library_kits: "30 × 25,000円 = 750,000円"
        AWS: "3,000円 × 30 = 90,000円"
        小計: "1,110,000円"

      MiSeq_外注:
        外注費: "30 × 150,000円 = 4,500,000円"
        送料: "30 × 2,000円 = 60,000円"
        小計: "4,560,000円"

      共通試薬: "570,000円"

      Phase_2合計: "6,240,000円"

    Phase_3（Month 11-12、実運用12サンプル）:

      MinION_内製:
        臨床サンプル12個:
          Flowcells: "1個 × 90,000円 = 90,000円"
          Library: "12 × 25,000円 = 300,000円"
          AWS: "12 × 3,000円 = 36,000円"
        小計: "426,000円"

      MiSeq_外注（QC用途のみ）:
        四半期1サンプル: "0サンプル（この期間）"
        コスト: "0円"

      Phase_3合計: "426,000円"

  #===== 人件費 =====
  Personnel_Costs:

    研究期間12ヶ月:

      Principal_Investigator:
        時間: "20%FTE"
        給与: "10,000,000円/年 × 0.2 = 2,000,000円"

      Bioinformatician:
        時間: "50%FTE"
        給与: "6,000,000円/年 × 0.5 = 3,000,000円"

      Lab_Technician（2名）:
        時間: "100%FTE"
        給与: "4,000,000円/年 × 2 = 8,000,000円"

      人件費合計: "13,000,000円"

  #===== その他費用 =====
  Other_Costs:

    AWS_Storage: "50,000円/年"

    消耗品: "500,000円"

    学会発表・論文: "1,000,000円"

    外注先契約管理: "100,000円"

    予備費（10%）: "1,400,000円"

  #===== 総予算（改訂版） =====
  Total_Budget:

    初期投資: "7,450,000円"

    運営費（12ヶ月）:
      試薬・外注費: "12,906,000円"
      人件費: "13,000,000円"
      AWS: "150,000円"
      その他: "1,600,000円"
      予備費: "1,400,000円"

    運営費合計: "29,056,000円"

    総予算: "36,506,000円"

  #===== コスト比較 =====
  Cost_Comparison:

    【改訂計画】MiSeq外注:
      初期投資: "7.45M円"
      運営費: "29.06M円"
      Total: "36.5M円"

    【従来計画】MiSeq購入:
      初期投資: "24.45M円"
      運営費: "26.01M円"
      Total: "50.5M円"

    削減額: "-14M円（28%削減）"

  #===== Phase 3以降のコスト（Year 2-） =====
  Post_Study_Cost:

    MinION_完全内製運用（24サンプル/年）:
      試薬: "25,000円 × 24 = 600,000円"
      Flowcell: "90,000円 × 2 = 180,000円"
      AWS: "3,000円 × 24 = 72,000円"
      人件費: "3,000,000円/年"
      小計: "3,852,000円/年"

    MiSeq_外注（四半期QC、4サンプル/年）:
      外注費: "150,000円 × 4 = 600,000円/年"

    Year_2以降年間総コスト: "4,452,000円/年"

    vs. MiSeq外注全量（24サンプル）:
      外注費: "150,000円 × 24 = 3,600,000円"
      前処理: "8,000円 × 24 = 192,000円"
      人件費: "500,000円"
      Total: "4,292,000円/年"

    差額: "+160,000円/年（MinIONやや高い）"

    但し付加価値:
      ✓ TAT: 48時間 vs 外注3週間
      ✓ 緊急対応能力
      ✓ PERV全ゲノム解析
      ✓ 完全内製（外注依存ゼロ）
      ✓ ノウハウ蓄積

  結論:
    "わずか+280,000円/年で完全内製化達成、
     初期投資は14M円削減。ROI極めて高い。"
```

----

## 9. 期待される成果（改訂版）

### 9.1 学術的成果

```yaml
Academic_Outcomes:

  #===== 論文発表 =====
  Publications:

    Target_1:
      Title: "Validation of Oxford Nanopore MinION for pathogen metagenomics in xenotransplantation donor pigs: Comparison with outsourced Illumina MiSeq"
      Journal: "Xenotransplantation (IF: 3.5)"
      Timeline: "Month 12（投稿）"
      強調点: "外注MiSeqをReference standardとした新規アプローチ"

    Target_2:
      Title: "Cost-effective in-house pathogen detection using MinION for clinical xenotransplantation"
      Journal: "Journal of Clinical Microbiology (IF: 5.0)"
      Timeline: "Month 15（投稿）"
      強調点: "内製化による経済性とTAT短縮"

  #===== 学会発表 =====
  Conference_Presentations:

    International_Xenotransplantation_Association:
      時期: "Year 2"
      形式: "Oral presentation"
      テーマ: "MinION内製化戦略"

    日本移植学会:
      時期: "Year 1-2"
      形式: "ポスター発表"

  #===== 知的財産 =====
  Intellectual_Property:

    可能性:
      - "MinION最適化プロトコル（ノウハウ）"
      - "AWS解析パイプライン（ソフトウェア）"
      - "外注併用バリデーション手法（方法論）"

    保護:
      - "営業秘密として管理"
      - "ソフトウェア著作権登録検討"
```

### 9.2 臨床的成果

```yaml
Clinical_Outcomes:

  #===== 完全内製化達成 =====
  Full_In_House_Capability:

    成果:
      ✓ "MinION単独検査体制確立（外注依存ゼロ）"
      ✓ "PMDA承認可能なバリデーションデータ"
      ✓ "SOP・解析パイプライン完備"
      ✓ "訓練されたスタッフ"
      ✓ "初期投資最小化（-17M円削減）"

    価値:
      - "外部委託完全脱却"
      - "迅速な意思決定（48時間TAT vs 外注3週間）"
      - "コスト完全管理"
      - "ノウハウ100%内部蓄積"
      - "秘密保持（外部流出リスクゼロ）"

  #===== 異種移植推進 =====
  Xenotransplantation_Advancement:

    貢献:
      1: "ドナーブタ安全性評価の迅速化"
      2: "移植前スクリーニングのTAT短縮"
      3: "PERV詳細解析の実現"
      4: "未知病原体検出能力"

    波及効果:
      - "日本の異種移植研究の自立性向上"
      - "PMDA審査の加速"
      - "患者への早期医療提供"

  #===== 経済的インパクト =====
  Economic_Impact:

    研究期間中の削減:
      初期投資削減: "-17M円（MiSeq購入不要）"
      研究総コスト削減: "-14M円（vs MiSeq購入計画）"

    運用期の経済性:
      Year_2以降:
        - MinION内製: "4.6M円/年"
        - MiSeq外注全量: "4.3M円/年"
        - 差額: "+0.3M円/年（MinIONやや高）"

      但し付加価値考慮:
        緊急対応価値: "移植遅延回避 +10M円/回相当"
        ノウハウ蓄積: "Priceless"

      結論: "経済的にほぼ互角、戦略的価値は内製が圧倒的優位"
```

----

## 10. 倫理的配慮・規制対応（改訂版）

### 10.1 倫理審査

```yaml
Ethical_Considerations:

  #===== 動物実験倫理 =====
  Animal_Ethics:

    該当性:
      - "ドナーブタからの採血"
      - "既存サンプル利用（新規動物実験なし）"

    必要手続き:
      - "動物実験委員会審査（施設依存）"
      - "3Rの原則遵守"

  #===== データ管理倫理 =====
  Data_Management_Ethics:

    外注先への提供:
      - "サンプルID匿名化"
      - "動物個体情報は非開示"
      - "研究目的のみ開示"

    秘密保持:
      - "外注先とNDA締結"
      - "データ破棄確認（研究終了後）"

  #===== 研究公正 =====
  Research_Integrity:

    データ改ざん防止:
      - "MinION生データ保存（ALCOA+）"
      - "外注先生データ受領・保存"
      - "解析パラメータ記録"

    選択的報告回避:
      - "研究計画事前登録（UMIN-CTR等）"
      - "全結果報告（陰性結果も）"

    利益相反:
      - "ONT社・外注先からの資金提供なし"
      - "独立研究として実施"
```

### 10.2 PMDA対応準備（改訂版）

```yaml
PMDA_Submission_Preparation:

  #===== 事前相談 =====
  Pre_Consultation:

    時期: "Month 2-3"

    目的:
      - "MinION内製 + MiSeq外注比較の妥当性確認"
      - "外注Reference standardの受入可否"
      - "必要バリデーション項目明確化"

    提出資料:
      - "研究計画書"
      - "MinION技術概要"
      - "外注先の信頼性証明（ISO認定等）"
      - "予備データ"

    想定される指摘:

      指摘_1: "外注Reference standardの信頼性は？"
      回答:
        - "ISO 15189認定施設"
        - "PMDA承認実績あり"
        - "生データ（FASTQ）受領し自社検証可能"

      指摘_2: "外注先でのプロトコル標準化は？"
      回答:
        - "事前に詳細プロトコル合意"
        - "各ラン QCレポート確認"
        - "逸脱時は再解析要求"

  #===== 最終提出資料（改訂版） =====
  Final_Submission_Package:

    時期: "Month 12"

    内容:

      1_バリデーション報告書:
        - "本研究の全結果"
        - "MinION vs MiSeq外注比較データ"
        - "統計解析詳細"
        - "結論"

      2_SOP:
        - "MinION内製検査法SOP"
        - "AWS解析パイプラインSOP"
        - "QC/QA手順"

      3_分析的性能:
        - "LOD"
        - "精度・再現性"
        - "PPA/NPA"
        - "MiSeq外注との一致率"

      4_参照資料:
        - "外注先情報（ISO認定証等）"
        - "外注プロトコル"
        - "文献リスト"

      5_実施体制:
        - "施設概要"
        - "MinION内製担当者資格"
        - "機器リスト"
```

----

## 11. 成功基準・中止基準（改訂版）

### 11.1 研究成功の定義

```yaml
Success_Criteria:

  #===== Phase 1成功基準 =====
  Phase_1_Success:

    必須（Must）:
      1: "3条件中、1条件以上でMiSeq外注との一致率 >90%達成"
      2: "最適条件のPPA >90%"
      3: "外注先の品質安定（QC合格率 >95%）"
      4: "重大な安全性問題なし"

    望ましい（Should）:
      1: "一致率 >95%達成"
      2: "コストがMiSeq外注の1/5以下"

    判定: "Must全達成でPhase 2移行"

  #===== Phase 2成功基準（Exit Criteria） =====
  Phase_2_Success_Exit_Criteria:

    Primary:
      1: "種同定一致率 >95%（95% CI下限 >90%）"
      2: "定量相関 R² >0.90"

    Secondary:
      3: "PPA >95%"
      4: "NPA >98%"
      5: "MinION内製 LOD <100 copies/mL"

    判定: "Primary全達成 + Secondary 3/5達成でMinION完全内製化承認"

  #===== 最終成功基準 =====
  Overall_Success:

    1: "Phase 2 Exit criteria達成"
    2: "MinION単独運用を2ヶ月実施（Phase 3）"
    3: "PMDA提出資料完成"
    4: "論文1報以上投稿"
    5: "外注費用を年間600,000円以下に削減（QC用途のみ）"

    判定: "全達成で研究大成功、完全内製化達成"
```

### 11.2 研究中止基準（改訂版）

```yaml
Stopping_Rules:

  #===== Phase 1中止基準 =====
  Phase_1_Stopping:

    Criterion_1_性能不足:
      条件: "30サンプル完了時、全条件で一致率 <85%"
      Action: "研究中止、MinION諦めてMiSeq外注継続"

    Criterion_2_外注先品質不良:
      条件: "QC不合格率 >20%、改善されず"
      Action:
        - "外注先変更"
        - "期間延長（+2ヶ月）"

    Criterion_3_予算超過:
      条件: "Phase 1で予算執行率 >120%"
      Action: "Phase 2サンプル数削減 or 研究縮小"

  #===== Phase 2中止基準 =====
  Phase_2_Stopping:

    Criterion_1_Exit_criteria未達:
      条件: "一致率 <90%"
      Action:
        Option_A: "Phase 1に戻り再最適化"
        Option_B: "MinION諦め、MiSeq外注継続"

    Criterion_2_外注費高騰:
      条件: "外注費が契約価格の150%超"
      Action:
        - "外注先交渉 or 変更"
        - "サンプル数削減"

  #===== 最終判断 =====
  Final_Decision:

    Phase_2失敗時の選択肢:
      Option_A: "MinION内製 + MiSeq外注ハイブリッド"
      Option_B: "MinION諦め、MiSeq外注単独"
      Option_C: "期間延長、追加研究"

    判断基準:
      - "達成度（90-95%なら継続価値）"
      - "予算余裕"
      - "PMDA見解"
      - "外注先継続可能性"
```

----

## 12. 結論とNext Steps（改訂版）

### 12.1 研究の意義

```
╔═══════════════════════════════════════════════════════════════╗
║  MinION検証研究の戦略的位置づけ（MiSeq外注版）                ║
╚═══════════════════════════════════════════════════════════════╝

【現状の課題】
  ❌ MinION: 新規技術、PMDA承認実績なし
  ❌ 単独使用の妥当性証明データ不足
  ✓ MiSeq: ゴールドスタンダードだが外注費高額（160,000円/サンプル）

【本研究の戦略（改訂版）】
  1. MiSeq外注をReference standardとして並行解析
  2. MinIONプロトコル・解析パイプライン最適化（Phase 1）
  3. 最適MinIONの独立検証（Phase 2）
  4. 統計的同等性証明 → MinION完全内製化（Phase 3）

【期待される成果】
  ✓ MinION単独でPMDA承認可能なバリデーションデータ取得
  ✓ 初期投資削減: -17M円（MiSeq装置不要）
  ✓ Phase 3以降、MiSeq外注ほぼゼロ（四半期1回QCのみ）
  ✓ TAT短縮: 48時間 vs MiSeq外注3週間
  ✓ PERV詳細解析能力（MinION長鎖readの優位性）
  ✓ 完全内製化達成（外部委託依存ゼロ）

【Exit Strategy】
  成功時: MinION完全内製運用、MiSeq外注は四半期QCのみ
  失敗時: MiSeq外注継続、MinIONは補助用途

【投資対効果（改訂版）】
  研究総予算: 36.5M円（12ヶ月）
    vs. MiSeq購入計画: 50.5M円
    削減額: -14M円（28%削減）

  初期投資: 7.5M円
    vs. MiSeq購入: 24.5M円
    削減額: -17M円（70%削減）

  Phase 3以降年間コスト: 4.6M円
    vs. MiSeq外注全量: 4.3M円
    差額: +0.3M円（ほぼ同等）

  付加価値:
    ✓ 完全内製化（外注依存ゼロ）
    ✓ TAT 48時間（vs 外注3週間）
    ✓ 緊急対応能力
    ✓ PERV全ゲノム解析
    ✓ ノウハウ100%内部蓄積
    ✓ 秘密保持（外部流出リスクゼロ）

  結論:
    "初期投資17M円削減し、運用コストほぼ同等で、
     完全内製化と迅速対応を実現。ROI極めて高い。"
```

### 12.2 即座のアクション

```yaml
Immediate_Actions:

  今週中:
    Priority_1:
      □ "本研究計画書（改訂版）を研究責任者・経営層に提出"
      □ "予算承認プロセス開始（36.5M円）"
      □ "MiSeq外注先RFP作成"

    Priority_2:
      □ "MinION調達状況確認"
      □ "AWS環境構築計画"

  今月中:
    □ "研究計画書最終化"
    □ "予算承認取得"
    □ "MiSeq外注先選定（RFP → 評価 → 契約）"
    □ "試薬発注（Phase 1分）"
    □ "研究チーム編成"

  Month_1開始:
    □ "外注先とプロトコル詳細合意"
    □ "パイロット試験（5サンプル）"
    □ "Phase 1サンプル収集開始"
```

----

**文書バージョン**: 2.0（MiSeq外注対応版）
**作成日**: 2025年10月8日
**改訂日**: 2025年10月8日
**次回更新**: Month 7（Phase 1完了時）、Month 10（Phase 2完了時）

**重要な変更点**:
- MiSeq装置購入を外注に変更
- 初期投資を17M円削減（7.5M円）
- 研究総予算を14M円削減（36.5M円）
- タイムラインを1ヶ月延長（外注TAT考慮）
- Phase 3以降、MiSeq外注を四半期1回QCのみに削減

**推奨アクション**: 本改訂版計画書を基に、MinION最適化研究を実施。12ヶ月後にMinION完全内製運用を実現し、外部委託依存から完全脱却。初期投資最小化とPMDA承認可能性を両立。
