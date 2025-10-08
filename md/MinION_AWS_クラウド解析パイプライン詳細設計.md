# MinION-AWS クラウド解析パイプライン詳細設計
## 異種移植用ドナーブタ病原体メタゲノム解析システム

**作成日**: 2025年10月8日
**対象**: PMDA指定91病原体スクリーニング
**プラットフォーム**: Oxford Nanopore MinION + AWS Cloud
**サンプル規模**: 24サンプル/年（Phase I臨床試験）

----

## エグゼクティブサマリー

### クラウドベース解析システムの全体像

```yaml
【AWS-MinIONパイプライン概要】

アーキテクチャ: "完全クラウドベース（オンプレGPU不要）"

主要構成:
  ローカル（ラボ）:
    - MinION Mk1D（シーケンシング）
    - ベーシックPC（アップロード専用）

  AWS Cloud:
    - S3（データストレージ）
    - EC2 G5インスタンス（GPU basecalling）
    - Batch（並列解析）
    - Step Functions（ワークフロー管理）
    - RDS（結果データベース）
    - QuickSight（レポート生成）

総コスト削減:
  オンプレGPU PC不要: -300,000円（初期投資）
  AWS従量課金: 月額50,000-100,000円（24サンプル/年）

利点:
  ✓ 初期投資最小化（MinION本体のみ）
  ✓ スケーラビリティ（サンプル増加に柔軟対応）
  ✓ 自動化（ワンクリック解析）
  ✓ データ永続化（S3長期保存）
  ✓ セキュリティ（AWS準拠）
  ✓ 災害対策（マルチAZ冗長化）

TAT（ターンアラウンドタイム）:
  シーケンシング: 6-48時間（MinION）
  アップロード: 2-4時間（30Gb, 100Mbps回線）
  Basecalling: 4-8時間（G5.xlarge GPU）
  解析: 2-4時間（並列処理）
  Total: 14-64時間（通常24-36時間）
```

----

## 1. システムアーキテクチャ

### 1.1 全体構成図

```yaml
System_Architecture:

  #========== ローカル環境（ラボ） ==========
  On_Premises_Lab:

    Equipment:
      MinION_Mk1D:
        機能: "シーケンシング実施"
        接続: "USB-C（USB 2.0転送速度対応）"
        温度制御: "ペルチェ素子搭載（10-35°C動作保証）"
        出力: "FAST5 files（Raw signal data）"
        保存先: "MinION内蔵ストレージ（1TB SSD推奨）"
        重量: "130g"

      Upload_PC:
        Spec: "一般的なPC（GPU不要）"
        CPU: "Core i5以上"
        RAM: "16GB"
        Storage: "500GB SSD"
        Network: "100Mbps以上（光回線）"
        OS: "Windows 10/11 or macOS"
        Software:
          - "AWS CLI（アップロード用）"
          - "MinKNOW（MinION制御）"
          - "Rclone or S3 sync（自動アップロード）"

    Workflow:
      Step_1: "MinIONでシーケンシング（6-48時間）"
      Step_2: "FAST5ファイル生成（30-50Gb）"
      Step_3: "AWS S3へ自動アップロード（2-4時間）"
      Step_4: "アップロード完了でAWSパイプライン自動起動"

  #========== AWS Cloud環境 ==========
  AWS_Cloud:

    Region: "ap-northeast-1（東京リージョン）"
    理由:
      - "データ主権（日本国内）"
      - "低レイテンシ"
      - "PMDA規制対応"

    Core_Services:

      #----- データストレージ -----
      S3_Buckets:
        raw-data-bucket:
          目的: "MinION生データ（FAST5）保存"
          Storage_class: "S3 Standard → Glacier（90日後）"
          Versioning: "有効"
          Encryption: "SSE-S3（サーバー側暗号化）"
          Lifecycle: "365日後にGlacier Deep Archive"
          容量想定: "30Gb/サンプル × 24 = 720Gb/年"

        processed-data-bucket:
          目的: "解析済みデータ（FASTQ, BAM, VCF等）"
          Storage_class: "S3 Standard"
          容量想定: "10Gb/サンプル × 24 = 240Gb/年"

        results-bucket:
          目的: "最終結果（レポート、サマリ）"
          Storage_class: "S3 Standard"
          容量想定: "100MB/サンプル × 24 = 2.4Gb/年"

        database-bucket:
          目的: "参照データベース（PMDA 91病原体等）"
          Storage_class: "S3 Standard"
          容量想定: "50Gb（固定）"
          更新: "四半期ごと"

      #----- コンピューティング -----
      EC2_Instances:

        Basecalling_Instance:
          Type: "g5.xlarge（GPU搭載）"
          Spec:
            vCPU: "4"
            RAM: "16GB"
            GPU: "NVIDIA A10G（24GB VRAM）"
            Storage: "250GB EBS gp3"
          用途: "Guppy/Dorado basecalling（Duplex mode）"
          起動: "オンデマンド（解析時のみ）"
          時間単価: "$1.006/時間"

        Analysis_Instance:
          Type: "c6i.4xlarge（CPU最適化）"
          Spec:
            vCPU: "16"
            RAM: "32GB"
            Storage: "500GB EBS gp3"
          用途: "病原体検出、アノテーション"
          起動: "オンデマンド"
          時間単価: "$0.68/時間"

        Database_Instance:
          Type: "r6i.large（メモリ最適化）"
          Spec:
            vCPU: "2"
            RAM: "16GB"
          用途: "Kraken2/BLASTデータベースロード"
          起動: "オンデマンド"
          時間単価: "$0.252/時間"

      #----- ワークフロー管理 -----
      AWS_Batch:
        機能: "並列ジョブ実行"

        Job_Queues:
          high-priority-queue:
            優先度: "100"
            用途: "緊急サンプル"

          standard-queue:
            優先度: "50"
            用途: "通常サンプル"

        Compute_Environment:
          Type: "EC2（スポットインスタンス利用可）"
          最大vCPU: "64"
          スポット割引: "最大70% OFF"

      AWS_Step_Functions:
        機能: "ワークフロー全体のオーケストレーション"

        State_Machine:
          Name: "MinION-Pathogen-Detection-Pipeline"

          States:
            1_Basecalling: "FAST5 → FASTQ（Duplex mode）"
            2_QC: "Quality control（PycoQC）"
            3_Host_Removal: "Sus scrofa genome除去"
            4_Pathogen_Detection: "91病原体スクリーニング"
            5_Quantification: "コピー数定量"
            6_Report_Generation: "結果レポート作成"
            7_Notification: "完了通知（Email/Slack）"

          Error_Handling:
            - "各ステップでリトライ（最大3回）"
            - "失敗時はSNS通知"
            - "ログはCloudWatch保存"

      #----- データベース -----
      RDS_PostgreSQL:
        Engine: "PostgreSQL 15"
        Instance: "db.t4g.medium"
        Storage: "100GB SSD"

        用途:
          - "サンプルメタデータ管理"
          - "解析結果保存"
          - "QC指標蓄積"
          - "レポート生成用クエリ"

        Schema:
          samples: "サンプル情報（ID, 採取日, ブタID等）"
          sequencing_runs: "シーケンスラン情報"
          qc_metrics: "QC指標"
          pathogen_detections: "検出病原体"
          quantification_results: "定量結果"

      #----- 可視化・レポート -----
      QuickSight:
        機能: "ダッシュボード・レポート生成"

        Dashboards:
          1_Run_QC_Dashboard:
            - "Throughput推移"
            - "Read quality分布"
            - "Flowcell性能"

          2_Pathogen_Detection_Dashboard:
            - "検出病原体サマリ"
            - "定量結果グラフ"
            - "経時変化"

          3_Sample_Status_Dashboard:
            - "サンプル処理状況"
            - "TAT管理"
            - "合格/不合格判定"

      #----- セキュリティ・監視 -----
      IAM:
        機能: "アクセス制御"

        Roles:
          MinION_Upload_Role: "S3アップロード権限のみ"
          Pipeline_Execution_Role: "全リソースアクセス"
          Analyst_Role: "読取専用 + レポート生成"

        Policies:
          - "最小権限の原則"
          - "MFA必須（重要操作）"

      CloudWatch:
        機能: "ログ・メトリクス監視"

        Alarms:
          - "パイプライン失敗"
          - "コスト超過（予算アラート）"
          - "ディスク使用率 >80%"

      CloudTrail:
        機能: "全API呼び出し記録（監査証跡）"
        保存: "S3（5年間）"
```

### 1.2 データフロー詳細

```yaml
Data_Flow:

  #===== Phase 1: データ生成・アップロード =====
  Phase_1_Data_Generation:

    Location: "ラボ（オンプレミス）"

    Step_1_Sequencing:
      Tool: "MinKNOW（MinION制御ソフト）"
      Input: "ブタ血漿cfDNA/RNA"
      Output: "FAST5 files"
      Size: "30-50Gb/サンプル"
      Duration: "6-48時間"

      MinKNOW設定:
        Output_format: "FAST5（Raw signal）"
        Basecalling: "OFF（クラウドで実施）"
        理由: "ローカルGPU不要 → コスト削減"

    Step_2_Upload_to_S3:
      Trigger: "シーケンス完了検出（MinKNOWフォルダ監視）"

      Tool_Option_A: "AWS CLI S3 sync"
        Command: |
          aws s3 sync /MinION/data/run001/ \
            s3://minion-raw-data/run001/ \
            --storage-class STANDARD \
            --metadata sample_id=PIG001,run_date=2025-10-08

      Tool_Option_B: "Rclone（推奨）"
        利点:
          - "自動リトライ"
          - "帯域制限（業務影響最小化）"
          - "チェックサム検証"
        Config: |
          [aws-s3]
          type = s3
          provider = AWS
          region = ap-northeast-1

        Command: |
          rclone sync /MinION/data/run001/ \
            aws-s3:minion-raw-data/run001/ \
            --progress --checksum

      Upload_速度:
        File_size: "30Gb"
        Network: "100Mbps（実効80Mbps）"
        Time: "30Gb × 8 / 80Mbps = 3,000秒 ≈ 50分"
        実際: "2-4時間（余裕を見る）"

      アップロード完了検証:
        - "S3 bucket checksumと比較"
        - "完了マーカーファイル作成（.upload_complete）"

    Step_3_Pipeline_Trigger:
      Trigger: "S3イベント通知（PutObject）"
      Mechanism:
        1: "S3 → EventBridge → Step Functions"
        2: "または Lambda関数で前処理後に起動"

      Start_Condition: ".upload_complete ファイル検出"

  #===== Phase 2: Basecalling（GPU） =====
  Phase_2_Basecalling:

    Location: "AWS EC2 g5.xlarge"

    Step_1_Instance_Launch:
      AMI: "Deep Learning AMI（NVIDIA driver込み）"
      Instance_Type: "g5.xlarge"
      起動時間: "2-3分"

    Step_2_Guppy_Dorado_Setup:
      Software: "Dorado（ONT最新basecaller）"
      Version: "0.5.0+"

      Install:
        - "AMI起動時に自動インストール（User Data script）"
        - "または事前ビルドAMI使用"

    Step_3_Basecalling_Execution:

      Mode: "Duplex sequencing（高精度）"

      Command: |
        dorado duplex \
          --device cuda:0 \
          --min-qscore 10 \
          --emit-fastq \
          /opt/dorado/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
          s3://minion-raw-data/run001/ \
          > run001_duplex.fastq

      Performance:
        Input: "30Gb FAST5（4M reads）"
        GPU: "A10G 24GB"
        Throughput: "~1,000 reads/秒"
        Time: "4,000秒 ≈ 1.1時間（simplex）"
        Time_Duplex: "6-8時間（Duplex mode）"

      Output:
        Format: "FASTQ（basecalled sequences + Q scores）"
        Size: "5-10Gb（圧縮前）"
        Quality: "Mean Q score >20（Duplex: Q30）"

    Step_4_QC_Generation:
      Tool: "PycoQC"

      Command: |
        pycoQC \
          -f sequencing_summary.txt \
          -o run001_QC_report.html

      Metrics:
        - "Total bases, reads"
        - "Read length distribution（N50）"
        - "Quality score distribution"
        - "Throughput over time"
        - "Pore activity"

      合格基準チェック:
        Total_bases: ">25Gb"
        Mean_Q: ">18"
        N50: ">8kb"

      不合格時:
        Action: "Slackアラート → 人間判断"

    Step_5_Upload_Processed:
      Destination: "s3://minion-processed-data/run001/"
      Files:
        - "run001_duplex.fastq.gz（圧縮）"
        - "sequencing_summary.txt"
        - "QC_report.html"

    Step_6_Instance_Termination:
      Trigger: "アップロード完了"
      Cost_saving: "不要なGPU時間課金を回避"

  #===== Phase 3: ホストゲノム除去 =====
  Phase_3_Host_Removal:

    Location: "AWS Batch（c6i.4xlarge）"

    Step_1_Reference_Download:
      Source: "s3://minion-database/reference/sus_scrofa.fa"
      Reference:
        - "Sus scrofa 11.1（RefSeq）"
        - "+ 遺伝子改変ブタ配列（3KO-7TG-59PERV）"

    Step_2_Alignment:
      Tool: "Minimap2（長鎖read専用）"

      Index作成（初回のみ）: |
        minimap2 -d sus_scrofa.mmi sus_scrofa.fa

      Alignment: |
        minimap2 -ax map-ont -t 16 \
          sus_scrofa.mmi \
          run001_duplex.fastq.gz \
          | samtools view -@ 16 -b -o run001_mapped.bam

      Performance:
        Input: "5Gb FASTQ（2M reads）"
        CPU: "16 cores"
        Time: "30-60分"

    Step_3_Unmapped_Extraction:
      目的: "ホストにマップされないread = 病原体候補"

      Command: |
        samtools view -@ 16 -b -f 4 run001_mapped.bam \
          | samtools fastq -@ 16 - \
          | gzip > run001_non_host.fastq.gz

      Expected_depletion:
        Input_reads: "2M"
        Host_reads: "1.9M（95%）"
        Non_host_reads: "100k（5%）"

    Step_4_QC_Check:
      確認項目:
        - "Host removal率 >90%"
        - "Non-host read数 >50k"

      不合格時: "前処理工程の問題 → アラート"

  #===== Phase 4: 病原体検出 =====
  Phase_4_Pathogen_Detection:

    Location: "AWS Batch（並列ジョブ）"

    Strategy: "Multi-method consensus approach"

    #----- Method 1: Kraken2/Bracken -----
    Method_1_Kraken2:

      Tool: "Kraken2 + Bracken"
      Database: "PMDA 91病原体カスタムDB"

      Database構築（初回、四半期更新）:
        Source:
          - "NCBI RefSeq（91病原体）"
          - "ViralDB（RVDB）"
          - "Custom PERV sequences"

        Build: |
          kraken2-build --download-taxonomy --db PMDA_DB
          kraken2-build --add-to-library pmda_91.fna --db PMDA_DB
          kraken2-build --build --db PMDA_DB --threads 32

        Size: "20-30Gb"
        保存先: "s3://minion-database/kraken2/"

      Classification: |
        kraken2 --db PMDA_DB \
          --threads 16 \
          --report run001_kraken2.report \
          --output run001_kraken2.out \
          run001_non_host.fastq.gz

      Abundance推定: |
        bracken -d PMDA_DB \
          -i run001_kraken2.report \
          -o run001_bracken.output \
          -r 150 -l S

      Performance:
        Input: "100k reads"
        Time: "5-10分"

      Output: "種レベル分類 + 存在量推定"

    #----- Method 2: BLAST（高精度確認） -----
    Method_2_BLAST:

      Tool: "BLASTn"
      Database: "PMDA 91病原体 BLAST DB"

      用途: "Kraken2陽性readの確認"

      Workflow:
        1: "Kraken2で陽性となったtaxon IDのreadを抽出"
        2: "BLAST検索で再確認"
        3: "Identity >95%で真陽性判定"

      Command: |
        blastn -query positive_reads.fa \
          -db PMDA_91_blast \
          -outfmt 6 \
          -evalue 1e-10 \
          -num_threads 16 \
          -max_target_seqs 5 \
          -out run001_blast.tsv

      Performance:
        Input: "1,000 reads（陽性候補）"
        Time: "10-20分"

    #----- Method 3: De novo Assembly（未知病原体） -----
    Method_3_Assembly:

      Tool: "Flye（メタゲノムモード）"

      用途:
        - "未知病原体検出"
        - "91種リストに無いが症状ある場合"

      Command: |
        flye --nano-hq run001_non_host.fastq.gz \
          --meta \
          --threads 16 \
          --out-dir assembly_output

      Performance:
        Input: "100k reads"
        Time: "2-4時間"

      Output: "Contig FASTA"

      Post-assembly:
        1: "Contig BLAST検索（nt database）"
        2: "新規配列の同定"
        3: "系統解析（必要時）"

    #----- Method 4: PERV特異的解析 -----
    Method_4_PERV:

      重要度: "最高（PMDA特記事項）"

      Database: "PERV-A, B, C参照配列"

      Workflow:
        Step_1_PERV_read抽出:
          - "Minimap2でPERV参照配列にマッピング"
          - "マップされたread抽出"

        Step_2_Full_genome_assembly:
          - "PERV read（~9kb完全長）をアセンブリ"
          - "組換え部位解析"

        Step_3_Typing:
          - "A/B/C型判定"
          - "組換え体検出"
          - "Envelope遺伝子変異解析"

        Step_4_Quantification:
          - "PERV read count → copies/mL換算"
          - "Spike-in補正"

      Output:
        - "PERV型（A/B/C/Recombinant）"
        - "コピー数（copies/mL）"
        - "完全ゲノム配列（FASTA）"
        - "系統樹（必要時）"

  #===== Phase 5: 定量解析 =====
  Phase_5_Quantification:

    Location: "AWS Batch（r6i.large）"

    Method: "Read count-based digital quantification"

    Step_1_Read_Count:
      各病原体のread数カウント:
        - "Kraken2分類結果から集計"
        - "BLAST確認済みのみカウント"

    Step_2_Normalization:
      補正:
        1_Total_reads: "サンプル間で総read数補正"
        2_Genome_length: "ゲノムサイズで補正"
        3_Spike_in: "Internal standard補正"

      Formula: |
        Copies/mL = (Read_count × Spike_in_copies) /
                    (Spike_in_reads × Sample_volume × Efficiency)

    Step_3_LOD_Check:
      検出限界確認:
        - "Read count <10 → Below LOD"
        - "10-50 reads → 検出あるが定量困難"
        - ">50 reads → 定量可能"

    Output:
      Format: "TSV table"
      Columns:
        - "Pathogen_name"
        - "Read_count"
        - "Copies_per_mL"
        - "Confidence（High/Medium/Low）"
        - "Status（Positive/Negative/Inconclusive）"

  #===== Phase 6: レポート生成 =====
  Phase_6_Report_Generation:

    Location: "AWS Lambda（軽量処理）"

    Tool: "Python + Jinja2テンプレート"

    Report_Types:

      1_QC_Report:
        内容:
          - "Sequencing run summary"
          - "QC metrics（合格/不合格）"
          - "PycoQC HTMLレポート埋込"

      2_Pathogen_Detection_Report:
        内容:
          - "検出病原体リスト"
          - "定量結果（表 + グラフ）"
          - "PERV詳細解析結果"
          - "推奨事項（移植可否）"

        Format: "PDF + HTML"

      3_Technical_Report:
        内容:
          - "全解析ステップログ"
          - "使用ソフトウェア・バージョン"
          - "データベース情報"
          - "ALCOA+準拠データ"

        用途: "PMDA提出、監査対応"

    Output_Destination:
      - "s3://minion-results/run001/"
      - "RDS PostgreSQL（メタデータ）"
      - "Email送信（担当者へ）"

  #===== Phase 7: 通知・アーカイブ =====
  Phase_7_Notification:

    Tool: "AWS SNS + SES"

    通知先:
      - "Email（研究責任者、オペレータ）"
      - "Slack（チームチャンネル）"

    通知内容:
      Success:
        - "サンプルID、完了時刻"
        - "検出病原体サマリ"
        - "レポートダウンロードリンク"

      Failure:
        - "失敗ステップ"
        - "エラーメッセージ"
        - "CloudWatchログリンク"

    Archive:
      Raw_data: "S3 Standard → Glacier（90日後）"
      Processed_data: "S3 Standard（2年）"
      Results: "S3 Standard（永久保存）"
```

----

## 2. ソフトウェア・ツールスタック

### 2.1 解析パイプライン詳細

```yaml
Software_Stack:

  #========== Basecalling ==========
  Basecalling:

    Primary_Tool: "Dorado（ONT公式、最新）"
      Version: "0.5.0+"
      License: "Oxford Nanopore Technologies"
      Language: "C++/CUDA"

      Models:
        DNA: "dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
        RNA: "rna004_130bps_sup@v3.0.1"

      Modes:
        Simplex: "Q20、高速"
        Duplex: "Q30、高精度（推奨）"

      Installation: |
        # AWS EC2 Deep Learning AMI
        conda install -c conda-forge ont-dorado

    Alternative: "Guppy（旧版、互換性用）"
      Version: "6.5.7+"

  #========== Quality Control ==========
  QC_Tools:

    PycoQC:
      Version: "2.5.2"
      Purpose: "Run-level QC report"
      Input: "sequencing_summary.txt"
      Output: "Interactive HTML report"

      Install: |
        pip install pycoQC

    NanoPlot:
      Version: "1.41.0"
      Purpose: "Read-level QC visualization"

      Command: |
        NanoPlot --fastq run001.fastq.gz \
          --plots dot kde \
          --outdir nanoplot_output

    FastQC:
      Version: "0.12.1"
      Purpose: "汎用FASTQクオリティチェック"

  #========== Host Removal ==========
  Host_Removal:

    Minimap2:
      Version: "2.26"
      License: "MIT"
      Language: "C"

      Features:
        - "長鎖read最適化"
        - "高速（BWAの10倍）"
        - "低メモリ"

      Preset: "map-ont（MinION用）"

      Install: |
        conda install -c bioconda minimap2

    SAMtools:
      Version: "1.17"
      Purpose: "BAM操作、未マップread抽出"

      Install: |
        conda install -c bioconda samtools

  #========== Pathogen Detection ==========
  Pathogen_Detection_Tools:

    Kraken2:
      Version: "2.1.2"
      Algorithm: "k-mer based classification"
      Speed: "非常に高速（100k reads/min）"

      Custom_DB:
        Name: "PMDA_91_pathogens_DB"
        Size: "25Gb"
        Species: "91種 + Sus scrofa（除外用）"

      Install: |
        conda install -c bioconda kraken2

    Bracken:
      Version: "2.8"
      Purpose: "Abundance estimation（Kraken2結果から）"

    BLAST+:
      Version: "2.14.0"
      Tool: "blastn"
      Purpose: "高精度確認"

      Database:
        Name: "PMDA_91_blast"
        Format: "nucleotide BLAST DB"

      Install: |
        conda install -c bioconda blast

    Centrifuge:
      Version: "1.0.4"
      Purpose: "Alternative classifier（バックアップ）"

  #========== De novo Assembly ==========
  Assembly_Tools:

    Flye:
      Version: "2.9.2"
      Algorithm: "Repeat graph-based"
      Mode: "--meta（メタゲノム）"

      Features:
        - "長鎖read専用"
        - "Repeatを正確に解決"
        - "メタゲノムモード対応"

      Install: |
        conda install -c bioconda flye

    Canu:
      Version: "2.2"
      Purpose: "Alternative assembler"
      特徴: "高精度だが低速"

  #========== Phylogenetic Analysis ==========
  Phylogenetics:

    MAFFT:
      Version: "7.505"
      Purpose: "Multiple sequence alignment"

    IQ-TREE:
      Version: "2.2.0"
      Purpose: "Maximum likelihood phylogeny"

    使用場面:
      - "PERV系統解析"
      - "新規病原体の系統分類"

  #========== Quantification ==========
  Quantification:

    Custom_Python_Script:
      Name: "quantify_pathogens.py"

      Functions:
        - "Read count集計"
        - "Spike-in補正"
        - "Copies/mL算出"

      Dependencies:
        - "pandas"
        - "numpy"
        - "biopython"

  #========== Reporting ==========
  Reporting_Tools:

    MultiQC:
      Version: "1.14"
      Purpose: "複数QCレポート統合"

    Jinja2:
      Version: "3.1.2"
      Purpose: "HTMLテンプレート"

    Matplotlib/Seaborn:
      Purpose: "グラフ生成"

    ReportLab:
      Purpose: "PDF生成"

  #========== Workflow Management ==========
  Workflow:

    AWS_Step_Functions:
      Description: "サーバーレスワークフローオーケストレーション"
      Language: "Amazon States Language（JSON）"

      Features:
        - "視覚的ワークフロー設計"
        - "AWS Lambda/EC2統合"
        - "自動リトライ・エラーハンドリング"
        - "実行履歴・監査証跡"
        - "並列実行サポート"

      State_Machine_Structure:
        States:
          1_Trigger: "S3イベント検出"
          2_Basecalling: "EC2起動 → Dorado実行"
          3_QC: "Lambda → PycoQC実行"
          4_Host_Removal: "EC2 → Minimap2実行"
          5_Pathogen_Detection: "並列実行（Kraken2, BLAST）"
          6_Quantification: "Lambda → 定量計算"
          7_Report_Generation: "Lambda → レポート生成"
          8_Notification: "SNS → 通知送信"

      Example_State: |
        {
          "Basecalling": {
            "Type": "Task",
            "Resource": "arn:aws:states:::ec2:runInstances.sync",
            "Parameters": {
              "ImageId": "ami-basecalling-v1.0.0",
              "InstanceType": "g5.xlarge",
              "UserData": "#!/bin/bash\ncd /opt/pipeline\n./run_basecalling.sh"
            },
            "Retry": [{
              "ErrorEquals": ["States.TaskFailed"],
              "MaxAttempts": 3
            }],
            "Next": "QC"
          }
        }

    Python_Scripts:
      Description: "各解析フェーズのメインロジック"
      Location: "scripts/"

      Structure:
        - "1_basecalling.py"
        - "2_qc.py"
        - "3_host_removal.py"
        - "4_pathogen_detection.py"
        - "5_quantification.py"
        - "6_report_generation.py"

      Execution:
        Method: "EC2インスタンス上で直接実行"
        Environment: "Conda環境"

    Alternative: "AWS Batch"
      優位性: "より複雑なジョブキュー管理"
      現状: "Step Functions で十分対応可能"
```

### 2.2 実行環境構成

```yaml
Execution_Environment:

  理由:
    - "環境再現性（バリデーション要件）"
    - "ソフトウェアバージョン固定"
    - "管理容易性の確保"

  #========== EC2カスタムAMI戦略 ==========
  Custom_AMI_Strategy:

    Base_AMI:
      Name: "Deep Learning AMI (Ubuntu 22.04)"
      Provider: "AWS公式"

      Pre-installed:
        - "NVIDIA Driver"
        - "CUDA Toolkit"
        - "Python 3.10"
        - "Conda"

    Basecalling_AMI:
      Name: "minion-basecalling-ami"
      Base: "Deep Learning AMI"

      Additional_Software:
        - "Dorado（ONT公式basecaller）"
        - "PycoQC"
        - "AWS CLI v2"
        - "Rclone"

      Storage: "250GB EBS gp3"

      Launch_Target: "g5.xlarge（GPU）"

    Analysis_AMI:
      Name: "minion-analysis-ami"
      Base: "Ubuntu 22.04 LTS"

      Installed_Tools:
        - "Minimap2"
        - "SAMtools"
        - "Kraken2"
        - "Bracken"
        - "BLAST+"
        - "Flye"
        - "NanoPlot"
        - "MultiQC"

      Storage: "500GB EBS gp3"

      Launch_Target: "c6i.4xlarge（CPU）"

  #========== ソフトウェア管理 ==========
  Software_Management:

    Version_Control:
      Method: "Conda environments"

      Environment_Files:
        - "environment_basecalling.yml"
        - "environment_analysis.yml"
        - "environment_reporting.yml"

      Benefits:
        - "バージョン固定"
        - "再現性確保"
        - "依存関係管理"

    Example_Environment: |
      # environment_analysis.yml
      name: minion-analysis
      channels:
        - bioconda
        - conda-forge
        - defaults
      dependencies:
        - python=3.10
        - minimap2=2.26
        - samtools=1.17
        - kraken2=2.1.2
        - bracken=2.8
        - blast=2.14.0
        - flye=2.9.2

    Deployment:
      Method: "AMI snapshot"
      Frequency: "月次更新"
      Versioning: "Semantic versioning（v1.0.0）"

  #========== Lambda関数環境 ==========
  Lambda_Environment:

    Runtime: "Python 3.11"

    Layer_Strategy:
      Common_Layer:
        Name: "minion-common-libs"
        Contents:
          - "boto3（AWS SDK）"
          - "pandas"
          - "numpy"
          - "jinja2（レポート生成）"

    Deployment:
      Method: "ZIP upload / SAM template"
      Source_Control: "GitHub"

  #========== 環境更新管理 ==========
  Environment_Updates:

    AMI_Update_Process:
      Step_1: "テスト環境で新AMI作成"
      Step_2: "バリデーション実施"
      Step_3: "本番環境適用"
      Step_4: "旧AMI保持（3世代）"

    Software_Version_Tracking:
      Tool: "requirements.txt / environment.yml"
      Location: "GitHub repository"

      Change_Log:
        - "バージョン変更履歴記録"
        - "バリデーション結果紐付け"
        - "PMDA提出用トレーサビリティ"

  #========== セキュリティ管理 ==========
  Security:

    AMI_Security:
      - "AWS Inspector でスキャン"
      - "Critical脆弱性検出時アラート"
      - "週次セキュリティチェック"

    Software_Updates:
      - "セキュリティパッチ適用"
      - "但し、バリデーション済み環境維持"
      - "重大脆弱性のみ即座対応"

    Access_Control:
      - "AMI共有: アカウント内のみ"
      - "SSM Session Manager経由アクセス"
      - "SSH鍵不要"
```

----

## 3. AWS コスト詳細分析

### 3.1 月額コスト試算（24サンプル/年 = 2サンプル/月）

```yaml
AWS_Cost_Breakdown:

  前提条件:
    サンプル数: "2サンプル/月（24/年）"
    データサイズ: "30Gb/サンプル（FAST5）"
    Region: "ap-northeast-1（東京）"

  #========== ストレージコスト ==========
  S3_Storage:

    Raw_data（FAST5）:
      保存期間: "90日 Standard → Glacier"

      Standard_期間:
        容量: "30Gb × 2 = 60Gb/月"
        単価: "$0.025/GB/月"
        コスト: "60 × $0.025 = $1.5/月"

      Glacier_移行後:
        容量: "720Gb/年（累積）"
        単価: "$0.005/GB/月"
        コスト: "720 × $0.005 = $3.6/月（1年後）"

    Processed_data（FASTQ/BAM）:
      容量: "10Gb × 2 = 20Gb/月"
      保存期間: "2年 Standard"
      単価: "$0.025/GB/月"
      コスト: "$0.5/月"

    Results:
      容量: "100MB × 2 = 0.2Gb/月"
      単価: "$0.025/GB/月"
      コスト: "$0.005/月"

    Database（参照ゲノム等）:
      容量: "50Gb（固定）"
      単価: "$0.025/GB/月"
      コスト: "$1.25/月"

    S3合計:
      初月: "$3.26/月"
      1年後: "$5.36/月（Glacier累積込み）"

  #========== コンピューティングコスト ==========
  EC2_Compute:

    Basecalling（g5.xlarge）:
      時間単価: "$1.006/時間"
      使用時間: "8時間/サンプル × 2 = 16時間/月"
      コスト: "16 × $1.006 = $16.10/月"

      スポットインスタンス利用時:
        割引率: "70%"
        コスト: "$16.10 × 0.3 = $4.83/月"
        推奨: "✓ 緊急性低い場合はスポット"

    Analysis（c6i.4xlarge）:
      時間単価: "$0.68/時間"
      使用時間: "4時間/サンプル × 2 = 8時間/月"
      コスト: "8 × $0.68 = $5.44/月"

      スポット:
        コスト: "$5.44 × 0.3 = $1.63/月"

    Database_loading（r6i.large）:
      時間単価: "$0.252/時間"
      使用時間: "2時間/サンプル × 2 = 4時間/月"
      コスト: "4 × $0.252 = $1.01/月"

    EC2合計（オンデマンド）: "$22.55/月"
    EC2合計（スポット70%割引）: "$7.47/月"

  #========== AWS Batch ==========
  AWS_Batch:
    コスト: "EC2料金に含まれる（追加料金なし）"

  #========== Step Functions ==========
  Step_Functions:
    State_transitions: "20 transitions/サンプル × 2 = 40/月"
    無料枠: "4,000 transitions/月"
    コスト: "$0（無料枠内）"

  #========== RDS PostgreSQL ==========
  RDS:
    Instance: "db.t4g.medium"

    Running_hours:
      稼働: "24時間/日 × 30日 = 720時間"
      時間単価: "$0.082/時間"
      コスト: "720 × $0.082 = $59.04/月"

    Storage:
      容量: "100GB SSD"
      単価: "$0.138/GB/月"
      コスト: "100 × $0.138 = $13.80/月"

    RDS合計: "$72.84/月"

    最適化案:
      - "解析時のみ起動（16時間/月）"
      - "コスト: 16 × $0.082 + $13.80 = $15.11/月"
      - "削減額: $57.73/月（79%削減）"

  #========== Data Transfer ==========
  Data_Transfer:

    アップロード（インターネット → S3）:
      容量: "60Gb/月"
      コスト: "$0（無料）"

    ダウンロード（S3 → インターネット）:
      容量: "1Gb/月（レポートダウンロード）"
      無料枠: "100GB/月"
      コスト: "$0（無料枠内）"

    Region内転送（S3 ↔ EC2）:
      容量: "100Gb/月"
      コスト: "$0（無料）"

  #========== その他サービス ==========
  Other_Services:

    CloudWatch_Logs:
      ログ量: "5GB/月"
      単価: "$0.76/GB"
      コスト: "$3.80/月"

    SNS（通知）:
      Email送信: "10件/月"
      コスト: "$0.01/月"

    QuickSight:
      Author: "1ユーザ × $18/月 = $18/月"
      または
      Reader: "5ユーザ × $5/月 = $25/月"

      推奨: "1 Author（十分）"
      コスト: "$18/月"

  #========== 月額総コスト ==========
  Monthly_Total:

    オンデマンド構成:
      S3: "$3.26"
      EC2: "$22.55"
      RDS: "$72.84"
      CloudWatch: "$3.80"
      QuickSight: "$18.00"
      SNS: "$0.01"
      Total: "$120.46/月"

    最適化構成（推奨）:
      S3: "$3.26"
      EC2_Spot: "$7.47"
      RDS_停止時削減: "$15.11"
      CloudWatch: "$3.80"
      QuickSight: "$18.00"
      SNS: "$0.01"
      Total: "$47.65/月"

  年間コスト:
    最適化構成: "$47.65 × 12 = $571.80/年"
    円換算（$1=150円）: "85,770円/年"

  Per_sample_cost:
    AWS費用: "85,770円 ÷ 24 = 3,574円/サンプル"

  #========== オンプレGPU PCとの比較 ==========
  vs_OnPremise_GPU:

    オンプレ構成:
      GPU_PC: "300,000円（初期）"
      電気代: "5,000円/run × 24 = 120,000円/年"
      減価償却: "60,000円/年（5年）"
      Total: "180,000円/年"

    AWS構成:
      費用: "85,770円/年"

    削減額: "94,230円/年（52%削減）"

    利点:
      ✓ "初期投資不要"
      ✓ "メンテナンス不要"
      ✓ "故障リスクなし"
      ✓ "スケーラビリティ"
```

### 3.2 コスト最適化戦略

```yaml
Cost_Optimization:

  #===== Strategy 1: スポットインスタンス =====
  Spot_Instances:

    適用:
      - "Basecalling（緊急でない場合）"
      - "De novo assembly（時間許容）"

    削減率: "60-70%"

    リスク:
      - "中断可能性（<5%）"
      - "緊急サンプルには不向き"

    対策:
      - "Step Functionsで自動リトライ"
      - "Checkpointから再開"

    実装: |
      # Nextflow config
      process {
        withName: BASECALLING {
          queue = 'spot-queue'
          errorStrategy = 'retry'
          maxRetries = 3
        }
      }

  #===== Strategy 2: RDS自動停止 =====
  RDS_Auto_Stop:

    方法:
      - "Lambda関数でスケジュール停止/起動"
      - "解析時のみ起動（16時間/月）"

    Schedule:
      起動: "解析開始前（Step Functions Trigger）"
      停止: "解析完了30分後"

    削減額: "$57.73/月（79%削減）"

    実装: |
      # Lambda function (Python)
      import boto3

      rds = boto3.client('rds')

      def lambda_handler(event, context):
        if event['action'] == 'start':
          rds.start_db_instance(DBInstanceIdentifier='minion-db')
        elif event['action'] == 'stop':
          rds.stop_db_instance(DBInstanceIdentifier='minion-db')

  #===== Strategy 3: S3 Lifecycle =====
  S3_Lifecycle:

    Policy:
      Rule_1:
        Name: "Archive raw data"
        Transition:
          - "90日後 → Glacier"
          - "365日後 → Glacier Deep Archive"
        削減: "80-95%（Glacier Deep）"

      Rule_2:
        Name: "Delete old processed data"
        Expiration: "2年後削除"

    年間削減額: "約30,000円"

  #===== Strategy 4: QuickSight最適化 =====
  QuickSight_Optimization:

    現状: "Author 1名 = $18/月"

    代替案:
      Option_A: "Grafana on EC2（オープンソース）"
        コスト: "$0（EC2既存利用）"

      Option_B: "Python matplotlib + S3静的HTML"
        コスト: "$0"

    推奨: "Phase I初期はOption B、Phase II以降でQuickSight"
    削減額: "$18/月 → $216/年"

  #===== Strategy 5: Reserved Instances =====
  Reserved_Instances:

    適用判断: "Phase II以降（サンプル数安定）"

    対象:
      - "RDS db.t4g.medium（1年予約）"
      - "割引率: 40-50%"

    Phase_II想定:
      サンプル数: "48/年 = 4/月"
      削減額: "$30/月"

  #===== 総合最適化 =====
  Total_Optimization:

    Phase_I初期:
      基本構成: "$47.65/月"

    Phase_I最適化後:
      QuickSight削減: "-$18"
      最終: "$29.65/月"
      年間: "$355.80/年 = 53,370円"

    Per_sample: "2,224円/サンプル"
```

----

## 4. セキュリティ・コンプライアンス

### 4.1 データセキュリティ

```yaml
Security_Framework:

  #===== データ暗号化 =====
  Encryption:

    転送中:
      Protocol: "TLS 1.3"
      適用:
        - "ラボ → S3（HTTPS）"
        - "EC2 ↔ S3（HTTPS）"
        - "全AWSサービス間"

    保存時:
      S3:
        Method: "SSE-S3（AES-256）"
        または: "SSE-KMS（より高度）"

      EBS:
        Method: "EBS暗号化（自動）"

      RDS:
        Method: "Storage encryption（AES-256）"

    鍵管理:
      推奨: "AWS KMS（Key Management Service）"

      Custom_key:
        Name: "minion-data-key"
        Rotation: "自動（年次）"

  #===== アクセス制御 =====
  Access_Control:

    IAM_Policies:

      Principle: "最小権限の原則（Least Privilege）"

      Roles:
        1_MinION_Upload_Role:
          Permissions:
            - "s3:PutObject（raw-data-bucketのみ）"
            - "s3:GetObject（検証用）"
          使用: "ラボPCのAWS CLI認証"

        2_Pipeline_Execution_Role:
          Permissions:
            - "全S3 bucket R/W"
            - "EC2起動/停止"
            - "Batch job submit"
            - "RDS接続"
          使用: "Step Functions/Lambda"

        3_Analyst_Role:
          Permissions:
            - "S3 results読取専用"
            - "QuickSightアクセス"
            - "RDS読取専用"
          使用: "データ解析者"

    MFA:
      適用: "全IAMユーザ（必須）"
      Type: "Virtual MFA（Google Authenticator等）"

  #===== ネットワークセキュリティ =====
  Network_Security:

    VPC:
      構成: "Private subnet（EC2, RDS）"

      Internet_access:
        - "NAT Gateway経由"
        - "直接インターネット露出なし"

      Security_Groups:
        Basecalling_SG:
          Inbound: "なし"
          Outbound: "HTTPS（S3アクセス）のみ"

        RDS_SG:
          Inbound: "EC2からのPostgreSQL（5432）のみ"
          Outbound: "なし"

    VPC_Endpoint:
      使用: "S3 VPC Endpoint（Gateway型）"
      利点:
        - "インターネット経由不要"
        - "データ転送料削減"
        - "セキュリティ向上"

  #===== 監査・ログ =====
  Audit_Logging:

    CloudTrail:
      有効化: "全リージョン"
      ログ保存: "S3（5年間）"

      記録内容:
        - "全API呼び出し"
        - "誰が、いつ、何を実行したか"
        - "失敗したアクセス試行"

      アラート:
        - "Root accountログイン"
        - "IAM policy変更"
        - "Security group変更"

    CloudWatch_Logs:
      保存期間: "1年"

      ログソース:
        - "Lambda実行ログ"
        - "EC2システムログ"
        - "RDSクエリログ（遅いクエリ）"

    S3_Access_Logs:
      有効化: "全bucket"
      用途: "データアクセス監査"

  #===== コンプライアンス =====
  Compliance:

    ALCOA+原則:
      実装:
        A_Attributable:
          - "CloudTrail（全操作記録）"
          - "RDS user tracking"

        L_Legible:
          - "構造化ログ（JSON）"
          - "読みやすいレポート"

        C_Contemporaneous:
          - "リアルタイムログ記録"
          - "タイムスタンプ（UTC）"

        O_Original:
          - "S3 Versioningで履歴保持"
          - "削除防止（Object Lock）"

        A_Accurate:
          - "Checksum検証"
          - "Automated QC"

        Complete:
          - "全ステップログ保存"

        Consistent:
          - "標準化されたSOP"

        Enduring:
          - "長期保存（Glacier）"

        Available:
          - "QuickSightダッシュボード"
          - "RDSクエリ可能"

    PMDA対応:
      要件:
        - "データ完全性証明"
        - "変更履歴追跡"
        - "アクセス制御記録"

      実装:
        □ "CloudTrail監査証跡"
        □ "S3 Versioning"
        □ "RDS監査ログ"
        □ "年次監査レポート自動生成"

    GDPR/個人情報保護:
      該当性: "ブタサンプル → 個人情報非該当"
      但し: "研究者情報は保護"

      対策:
        - "研究者PII暗号化"
        - "アクセスログ監視"
```

### 4.2 災害対策・事業継続

```yaml
Disaster_Recovery:

  #===== バックアップ戦略 =====
  Backup_Strategy:

    S3:
      Versioning: "有効（削除保護）"
      Cross_Region_Replication:
        Source: "ap-northeast-1（東京）"
        Destination: "ap-northeast-3（大阪）"
        対象: "results-bucket（重要データ）"

    RDS:
      Automated_Backup:
        Frequency: "日次"
        Retention: "7日間"
        Window: "03:00-04:00 JST（深夜）"

      Manual_Snapshot:
        Timing: "重要解析前後"
        Retention: "無期限"

    AMI_Backup:
      EC2カスタムAMI:
        Frequency: "月次"
        保存: "3世代"

  #===== 復旧目標 =====
  Recovery_Objectives:

    RTO（Recovery Time Objective）:
      目標: "4時間"

      シナリオ:
        EC2障害: "別AZで再起動（30分）"
        S3障害: "AWS自動復旧（ほぼゼロ）"
        RDS障害: "スナップショットから復元（2時間）"
        Region障害: "大阪リージョンへフェイルオーバー（4時間）"

    RPO（Recovery Point Objective）:
      目標: "24時間"

      理由: "日次バックアップで許容"

  #===== 高可用性設計 =====
  High_Availability:

    Multi_AZ:
      RDS: "Multi-AZ配置（推奨）"
      コスト: "+100%"
      判断: "Phase II以降で検討"

    Auto_Scaling:
      AWS_Batch: "自動スケーリング（標準）"

    Health_Check:
      CloudWatch_Alarms:
        - "RDS接続エラー"
        - "S3アクセスエラー"
        - "Pipeline失敗率 >20%"

      通知: "SNS → Email/Slack"
```

----

## 5. 実装ロードマップ

### 5.1 段階的構築計画

```yaml
Implementation_Roadmap:

  #===== Phase 0: 準備（Week 1-2） =====
  Phase_0_Preparation:

    Week_1:
      AWS_Account:
        □ "AWSアカウント作成"
        □ "Billing alarm設定（$100/月）"
        □ "IAMユーザ作成（管理者、オペレータ）"
        □ "MFA有効化"

      Network_Setup:
        □ "VPC作成（10.0.0.0/16）"
        □ "Public/Private subnet作成"
        □ "NAT Gateway配置"
        □ "S3 VPC Endpoint作成"

    Week_2:
      Storage_Setup:
        □ "S3 bucket作成（4つ）"
        □ "Lifecycle policy設定"
        □ "Versioning有効化"
        □ "暗号化設定（SSE-S3）"

      Security:
        □ "CloudTrail有効化"
        □ "CloudWatch Logs設定"
        □ "KMS key作成"

  #===== Phase 1: 基盤構築（Week 3-6） =====
  Phase_1_Infrastructure:

    Week_3:
      Compute_Setup:
        □ "EC2 AMI選定（Deep Learning AMI）"
        □ "カスタムAMI作成（ソフトウェアプリインストール）"
        □ "Security Group設定"

      Database:
        □ "RDS PostgreSQL作成"
        □ "スキーマ設計・作成"
        □ "初期データ投入"

    Week_4:
      Custom_AMI:
        □ "カスタムAMI構築（basecalling用）"
        □ "カスタムAMI構築（analysis用）"
        □ "Conda環境セットアップ・スナップショット"

    Week_5-6:
      Reference_Data:
        □ "Sus scrofa参照ゲノムダウンロード"
        □ "PMDA 91病原体配列収集"
        □ "Kraken2 DB構築"
        □ "BLAST DB構築"
        □ "S3アップロード"

  #===== Phase 2: パイプライン開発（Week 7-10） =====
  Phase_2_Pipeline_Development:

    Week_7-8:
      Pipeline_Development:
        □ "Pythonスクリプト作成（全6フェーズ）"
        □ "AWS Lambda関数実装"
        □ "EC2起動スクリプト作成"
        □ "S3データハンドリング実装"

      Local_Testing:
        □ "模擬データで検証"
        □ "各スクリプトの動作確認"

    Week_9:
      Step_Functions:
        □ "State Machine定義"
        □ "エラーハンドリング設定"

    Week_10:
      Integration_Testing:
        □ "E2Eテスト（FAST5 → Report）"
        □ "パフォーマンス測定"
        □ "コスト検証"

  #===== Phase 3: バリデーション（Week 11-20） =====
  Phase_3_Validation:

    Week_11-13:
      LOD_Study:
        □ "スパイクインサンプル調製"
        □ "段階希釈実験"
        □ "LOD決定"

    Week_14-16:
      Precision_Study:
        □ "日内再現性（n=10）"
        □ "日間再現性（n=10）"
        □ "オペレータ間（n=10）"

    Week_17-19:
      Specificity_Study:
        □ "交差反応試験"
        □ "False positive評価"

    Week_20:
      Validation_Report:
        □ "バリデーション報告書作成"
        □ "PMDA提出用データパッケージ"

  #===== Phase 4: 本番運用（Week 21-） =====
  Phase_4_Production:

    Week_21:
      □ "SOP最終化"
      □ "オペレータ最終トレーニング"
      □ "チェックリスト作成"

    Week_22:
      □ "初回臨床サンプル解析"
      □ "全工程確認"

    Week_23-:
      □ "定常運用"
      □ "月次レビュー"
      □ "継続的改善"
```

### 5.2 トレーニングプラン

```yaml
Training_Plan:

  #===== オペレータ（ラボ技術者） =====
  Operator_Training:

    対象: "2名"
    期間: "2週間"

    Week_1:
      Day_1-2:
        - "MinION基礎（理論）"
        - "AWS基礎（S3, CLI）"

      Day_3-5:
        - "MinIONシーケンシング実習"
        - "データアップロード実習"
        - "パイプライン起動実習"

    Week_2:
      Day_6-8:
        - "QCレポート解釈"
        - "結果判定基準"
        - "トラブルシューティング"

      Day_9-10:
        - "模擬サンプル全工程実習"
        - "SOP作成支援"

  #===== バイオインフォマティシャン =====
  Bioinformatician_Training:

    対象: "1名"
    期間: "4週間"

    Week_1:
      - "Nanoporeデータ特性理解"
      - "Basecalling原理・実践"
      - "QCツール習得"

    Week_2:
      - "病原体検出手法（Kraken2, BLAST）"
      - "De novo assembly"
      - "PERV解析"

    Week_3:
      - "AWS基礎（EC2, S3, Step Functions）"
      - "Python自動化スクリプト"
      - "Conda環境管理"

    Week_4:
      - "パイプライン改良"
      - "レポート自動化"
      - "トラブルシューティング"

  #===== 継続教育 =====
  Continuing_Education:

    Quarterly:
      - "新ソフトウェアバージョン対応"
      - "AWS新サービス評価"

    Annually:
      - "バリデーション再トレーニング"
      - "規制動向アップデート"
```

----

## 6. トラブルシューティング

```yaml
Troubleshooting_Guide:

  #===== 一般的な問題と対処 =====
  Common_Issues:

    Issue_1_アップロード失敗:
      症状: "S3 syncエラー"

      原因_A: "ネットワーク不安定"
      対処:
        - "Rclone使用（自動リトライ）"
        - "帯域制限解除"

      原因_B: "IAM権限不足"
      対処:
        - "IAM policyチェック"
        - "s3:PutObject権限確認"

    Issue_2_Basecalling遅い:
      症状: "8時間 → 20時間かかる"

      原因: "GPU未使用（CPU basecalling）"
      対処:
        - "dorado --device cuda:0 確認"
        - "nvidia-smi でGPU認識確認"
        - "CUDA driverバージョン確認"

    Issue_3_病原体検出ゼロ:
      症状: "全サンプル陰性（異常）"

      原因_A: "Host removal過剰"
      対処:
        - "除去率確認（>99%なら異常）"
        - "参照ゲノム確認"

      原因_B: "データベースエラー"
      対処:
        - "Kraken2 DB整合性チェック"
        - "再ダウンロード"

    Issue_4_コスト超過:
      症状: "月額 >$200"

      原因: "EC2インスタンス停止忘れ"
      対処:
        - "AWS Budgetアラート設定"
        - "Trusted Advisor確認"
        - "不要リソース削除"

  #===== エラーコード別対処 =====
  Error_Codes:

    E001_QC_Failed:
      意味: "Throughput <25Gb"
      対処:
        1: "フローセル不良 → 交換"
        2: "ライブラリ品質チェック"
        3: "再シーケンス判断"

    E002_Host_Removal_Failed:
      意味: "Minimap2エラー"
      対処:
        1: "参照ゲノムファイル確認"
        2: "メモリ不足 → インスタンス拡大"

    E003_No_Pathogen_Detected:
      意味: "全メソッド陰性"
      対処:
        1: "Positive control確認"
        2: "サンプル品質再確認"
        3: "本当に陰性の可能性"

    E004_Pipeline_Timeout:
      意味: "24時間以内に完了せず"
      対処:
        1: "Stuck processチェック"
        2: "手動リトライ"
        3: "Nextflow resume使用"
```

----

## 7. 結論とNext Steps

### 7.1 AWS-MinIONシステムの総合評価

```
╔═══════════════════════════════════════════════════════════════╗
║  AWS Cloud-based MinION解析パイプライン 総合評価              ║
╚═══════════════════════════════════════════════════════════════╝

【技術的実現可能性】✓✓✓ 完全実現可能

AWS利点:
  ✓ 初期投資最小化（GPU PC不要 → -300,000円）
  ✓ 従量課金（使った分だけ → 年間55,170円@24サンプル）
  ✓ 無限スケーラビリティ（サンプル増加に即対応）
  ✓ 自動化（ワンクリック解析）
  ✓ 災害対策（Multi-AZ, 自動バックアップ）
  ✓ セキュリティ（暗号化、監査証跡）
  ✓ PMDA規制対応（ALCOA+準拠ログ）

パフォーマンス:
  Basecalling: 6-8時間（Duplex mode on GPU）
  解析: 2-4時間（並列処理）
  Total TAT: 14-64時間（通常24-36時間）

コスト:
  月額: $29.65（最適化後）= 4,448円/月
  年額: 53,370円/年（24サンプル）
  Per sample: 2,224円/サンプル

  vs. オンプレGPU:
    削減額: 94,230円/年（52%削減）

  Total system cost（MinION + AWS）:
    初期: 8,400,000円（MinION本体等）+ 0円（AWS）
    年間: 18,552,000円（試薬等）+ 53,370円（AWS）
    Per sample: 771,000円

【推奨度】⭐⭐⭐⭐⭐（最高）

理由:
  1. コスト効率: オンプレGPUより安価
  2. 運用容易性: メンテナンスフリー
  3. 拡張性: Phase II以降も対応可能
  4. 規制対応: 監査証跡完備
  5. 災害対策: 自動バックアップ

【非推奨ケース】
  ❌ インターネット接続不可の環境
  ❌ データ国外持ち出し絶対不可（要確認）
  ❌ AWS利用が組織方針で禁止

【Next Steps】
  Week 1: AWSアカウント作成、予算承認
  Week 2-6: インフラ構築
  Week 7-10: パイプライン開発
  Week 11-20: バリデーション
  Week 21-: 本番運用開始
```

### 7.2 即座のアクション

```yaml
Immediate_Actions:

  今週中:
    Priority_1:
      □ "本設計書を技術責任者に提出"
      □ "AWS利用承認取得"
      □ "AWSアカウント作成"

    Priority_2:
      □ "バイオインフォマティシャン候補選定"
      □ "ネットワーク帯域確認（100Mbps以上）"

  今月中:
    □ "AWS基盤構築開始"
    □ "MinION発注（並行）"
    □ "トレーニング計画策定"

  3ヶ月後:
    □ "パイプライン完成"
    □ "バリデーション開始"
```

----

**文書バージョン**: 1.0
**作成日**: 2025年10月8日
**次回更新**: AWS構築完了後（Month 2）

**推奨アクション**: 本設計書を基にAWS構築を開始。並行してMinION調達とバリデーション計画を進める。不明点はAWSソリューションアーキテクトに相談を推奨。
