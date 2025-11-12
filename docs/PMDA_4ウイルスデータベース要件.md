# PMDA 4ウイルスデータベース要件

**バージョン:** 1.0
**作成日:** 2025年11月12日
**目的:** PMDA指定4ウイルスの高感度検出のための参照ゲノムデータベース要件

---- 

## 概要

本文書では、以下の検出に必要な参照ゲノムデータベースを規定します：
1. ポリオーマウイルス（dsDNA）
2. ハンタウイルス（ssRNA-、3分節）
3. 東部ウマ脳炎ウイルス - EEEV（ssRNA+）
4. 豚スプーマウイルス（レトロウイルス、プロウイルスDNA）

すべてのデータベースは `/mnt/efs/databases/pmda/2024.1/` に保存されます。

---- 

## 1. ポリオーマウイルスデータベース

### 1.1 ディレクトリ構造
```
/mnt/efs/databases/pmda/2024.1/polyomavirus/
├── polyoma_all.fasta           # 全ポリオーマウイルスゲノム
├── polyoma_all.mmi             # Minimap2インデックス
├── polyoma_kraken2/            # Kraken2データベース
│   ├── hash.k2d
│   ├── opts.k2d
│   └── taxo.k2d
├── polyoma_blast/              # BLASTデータベース
│   ├── polyoma.nhr
│   ├── polyoma.nin
│   └── polyoma.nsq
└── metadata.tsv                # 種、アクセッション番号、系統
```

### 1.2 必要な参照ゲノム

| 種                    | アクセッション番号  | 長さ（bp） | 備考                 |
| -------------------- | ---------- | ------ | ------------------ |
| BKポリオーマウイルス          | NC\_001538 | 5,153  | ヒトポリオーマウイルス、リファレンス |
| JCポリオーマウイルス          | NC\_001699 | 5,130  | ヒトポリオーマウイルス        |
| SV40                 | NC\_001669 | 5,243  | サルウイルス40           |
| Sus scrofaポリオーマウイルス2 | MH381769   | 5,185  | **重要：豚特異的**        |
| メルケル細胞ポリオーマウイルス      | NC\_010277 | 5,387  | クロスリファレンス          |
| カロリンスカ研究所ポリオーマウイルス   | NC\_018102 | 5,024  | クロスリファレンス          |

### 1.3 データベース構築

```bash
# NCBIから配列をダウンロード
esearch -db nucleotide -query "Polyomaviridae[Organism] AND complete genome" | \
  efetch -format fasta > polyoma_all.fasta

# Minimap2インデックスを構築
minimap2 -d polyoma_all.mmi polyoma_all.fasta

# Kraken2データベースを構築
kraken2-build --download-taxonomy --db polyoma_kraken2/
kraken2-build --add-to-library polyoma_all.fasta --db polyoma_kraken2/
kraken2-build --build --db polyoma_kraken2/

# BLASTデータベースを構築
makeblastdb -in polyoma_all.fasta -dbtype nucl -out polyoma_blast/polyoma
```

### 1.4 検出パラメータ

- **マッピングツール:** Minimap2（`-ax map-ont`）
- **最小リード長:** 100 bp
- **マッピング品質閾値:** Q10
- **信頼度閾値:** ≥100リード かつ ≥10×カバレッジ
- **交差反応性チェック:** PERV、豚パルボウイルスとの区別

---- 

## 2. ハンタウイルスデータベース

### 2.1 ディレクトリ構造
```
/mnt/efs/databases/pmda/2024.1/hantavirus/
├── hantavirus_L_segment.fasta
├── hantavirus_M_segment.fasta
├── hantavirus_S_segment.fasta
├── hantavirus_all.fasta        # 連結L+M+S
├── hantavirus_all.mmi
├── hantavirus_kraken2/
├── amplicon_primers/
│   ├── hantavirus_pool_A.bed   # プライマープールA座標
│   ├── hantavirus_pool_B.bed   # プライマープールB座標
│   └── primer_sequences.tsv
└── metadata.tsv
```

### 2.2 必要な参照ゲノム

| 種                   | セグメント | アクセッション番号  | 長さ（nt） | 地理的地域         |
| ------------------- | ----- | ---------- | ------ | ------------- |
| **ハンタウイルス**         | L     | NC\_005222 | 6,533  | アジア（韓国、中国、日本） |
|                     | M     | NC\_005219 | 3,651  |               |
|                     | S     | NC\_005218 | 1,696  |               |
| **ソウルウイルス**         | L     | NC\_005238 | 6,544  | 世界的（都市ネズミ）    |
|                     | M     | NC\_005236 | 3,651  |               |
|                     | S     | NC\_005237 | 1,715  |               |
| **ドブラバ-ベオグラードウイルス** | L     | NC\_005233 | 6,550  | ヨーロッパ         |
|                     | M     | NC\_005234 | 3,643  |               |
|                     | S     | NC\_005235 | 1,700  |               |
| **プーマラウイルス**        | L     | NC\_005224 | 6,530  | ヨーロッパ         |
|                     | M     | NC\_005223 | 3,623  |               |
|                     | S     | NC\_005222 | 1,764  |               |

### 2.3 増幅プライマー設計（高感度用）

**基準:** Kim et al., Viruses 2021（ハンタウイルスマルチプレックスPCR）

```
合計アンプリコン数：36（L: 20、M: 11、S: 5）
アンプリコンサイズ：400 bp
オーバーラップ：100 bp
プライマー設計ツール：Primal Scheme（ARTIC Network）
```

**プライマープールAの例（18アンプリコン）:**
```tsv
名前	配列	長さ	Tm	セグメント	開始	終了
HTNV_L_1_LEFT	AGCTGACACATCAGTGTGCC	20	60.2	L	1	400
HTNV_L_1_RIGHT	GGTCAACAGTACCTGCCTAG	20	60.1	L	380	780
HTNV_L_3_LEFT	TGGACCTGATACCAAGGTCC	20	60.3	L	701	1100
...
```

### 2.4 検出パラメータ

- **3セグメント一致必須:** L かつ M かつ S すべて検出
- **セグメント当たり最小リード数:** 50リード
- **マッピングツール:** Minimap2またはKraken2
- **アンプリコンマッピング:** 期待されるアンプリコン座標にアライン
- **分節ゲノムアセンブリ:** medakaまたはNanopolishでコンセンサス作成

---- 

## 3. 東部ウマ脳炎ウイルス（EEEV）/アルファウイルスデータベース

### 3.1 ディレクトリ構造
```
/mnt/efs/databases/pmda/2024.1/alphavirus/
├── alphavirus_all.fasta        # PMDA指定全アルファウイルス
├── alphavirus_all.mmi
├── eeev_references.fasta       # EEEV系統解析用
├── eeev_phylogeny/
│   ├── eeev_alignment.fasta
│   └── eeev_tree.nwk
├── alphavirus_kraken2/
└── metadata.tsv
```

### 3.2 必要な参照ゲノム

**PMDA指定アルファウイルス（厚労省異種移植指針）:**

| ウイルス                    | PMDA行  | アクセッション番号  | 長さ（nt） | Poly(A)+ |
| ----------------------- | ------ | ---------- | ------ | -------- |
| **東部ウマ脳炎ウイルス（EEEV）**    | 32/55  | NC\_003899 | 11,841 | あり       |
| **西部ウマ脳炎ウイルス（WEEV）**    | 32/55  | NC\_003908 | 11,722 | あり       |
| **ベネズエラウマ脳炎ウイルス（VEEV）** | 33/56  | NC\_001449 | 11,441 | あり       |
| **ゲタウイルス**              | 18/41  | NC\_003696 | 11,680 | あり       |
| チクングニアウイルス              | リファレンス | NC\_004162 | 11,805 | あり       |
| ロスリバーウイルス               | リファレンス | NC\_001544 | 11,659 | あり       |

### 3.3 EEEV系統リファレンス

系統解析による系統割り当て（北米 vs 南米）のため：

| 系統  | 代表株              | アクセッション番号  | 地理的起源    |
| --- | ---------------- | ---------- | -------- |
| 北米  | EEEV FL93-939    | NC\_003899 | 米国（フロリダ） |
| 南米  | EEEV BeAn 436284 | KP765787   | ブラジル     |
| 南米  | EEEV Argentina   | FJ402886   | アルゼンチン   |

### 3.4 検出パラメータ

- **マッピングツール:** Minimap2
- **最小リード数:** ≥100リード
- **カバレッジ:** ≥10×平均深度
- **種割り当て:** ベストヒットリファレンス
- **系統確認:** nsP1またはE1遺伝子をアライン → FastTree
- **交差反応性:** EEEV、WEEV、VEEV、ゲタウイルス間の区別

---- 

## 4. 豚スプーマウイルスデータベース

### 4.1 ディレクトリ構造
```
/mnt/efs/databases/pmda/2024.1/spumavirus/
├── sfv_pol_gene.fasta          # サルスプーマウイルスpol遺伝子
├── ffv_pol_gene.fasta          # ネコスプーマウイルスpol遺伝子
├── bfv_pol_gene.fasta          # ウシスプーマウイルスpol遺伝子
├── spumavirus_all_pol.fasta    # 全スプーマウイルスpol遺伝子
├── spumavirus_pol.mmi
├── perv_pol_gene.fasta         # 識別用
├── spumavirus_blast/           # 配列確認用BLAST
└── metadata.tsv
```

### 4.2 必要な参照配列

**重要:** NCBIに豚スプーマウイルス参照ゲノムは存在しません。異種間リファレンスを使用します。

| 種               | 遺伝子 | アクセッション番号                 | 長さ（bp）  | 豚との相同性（予測） |
| --------------- | --- | ------------------------- | ------- | ---------- |
| サルスプーマウイルス（SFV） | pol | NC\_001364 (nt 2000-5000) | \~3,000 | 30-50%（縮重） |
| ネコスプーマウイルス（FFV） | pol | NC\_001871 (nt 2100-5100) | \~3,000 | 30-50%     |
| ウシスプーマウイルス（BFV） | pol | NC\_001831 (nt 2050-5050) | \~3,000 | 30-50%     |
| **PERV（識別用）**   | pol | AF038600 (nt 4500-6500)   | \~2,000 | **豚内因性**   |

### 4.3 ネステッドPCRプライマー配列

**スプーマウイルスpol遺伝子（逆転写酵素ドメイン）を標的とする縮重プライマー:**

```
# アウタープライマー（第1 PCR）
FV-pol-F1: 5'-GGNCARATHGGNATGTTYGG-3' (縮重度: 96倍)
FV-pol-R1: 5'-CCRTCNCCRAANCCRTC-3' (縮重度: 64倍)
期待産物: ~800 bp

# インナープライマー（第2 PCR）
FV-pol-F2: 5'-ATHGGNCARGGNTTYACNAC-3' (縮重度: 96倍)
FV-pol-R2: 5'-GTRTCNGTYTTRTCNCC-3' (縮重度: 64倍)
期待産物: ~400 bp

# 縮重コード
N = A/T/G/C (4倍)
R = A/G (2倍、プリン)
Y = C/T (2倍、ピリミジン)
H = A/T/C (3倍、Gを除く)
```

### 4.4 検出および確認戦略

1. **ネステッドPCR**（PBMCゲノムDNAから）
   2. 標的：pol遺伝子（最も保存された領域）
   3. LOD：1-10 copies/10⁵ PBMC

2. **アガロースゲル確認**
   2. 期待されるバンド：\~400 bp
   3. 陰性コントロール：Spumavirus陰性PBMC DNA
   4. 陽性コントロール：SFV陽性霊長類DNA（利用可能な場合）

3. **Sangerシーケンシング（必須）**
   2. 両鎖をシーケンス
   3. NCBI ntデータベースに対してBLAST
   4. 同一性確認：Spumavirusと\>70%、PERVと\<50%

4. **系統解析**
   2. SFV/FFV/BFV/PERV pol遺伝子とアライン
   3. Spumaretrovirinaeとのクラスタリング確認（Gammaretrovirus/PERVではない）

---- 

## 5. 品質管理と検証

### 5.1 データベース検証チェックリスト

- [ ]() すべての参照配列がダウンロード・検証済（MD5チェックサム）
- [ ]() 重複配列なし（アクセッション番号で重複除去）
- [ ]() 分類ID正しく割り当て（NCBI Taxonomy）
- [ ]() Minimap2インデックス構築・テスト済
- [ ]() Kraken2データベース構築・テスト済
- [ ]() BLASTデータベースフォーマット・テスト済
- [ ]() メタデータファイル完全（種、アクセッション番号、長さ、備考）

### 5.2 検出検証（In Silico）

**陽性コントロール（合成リード）:**
- 模擬MinIONリードを生成（BadreadまたはNanoSim）
- 陰性豚血漿メタゲノムに添加
- 濃度：10、50、100、500、1,000 copies/mL
- ≥50 copies/mLでの検出を検証

**陰性コントロール:**
- 病原体フリー豚血漿メタゲノム
- 偽陽性検出なし

**交差反応性テスト:**
- ポリオーマウイルス vs 豚パルボウイルス、PERV
- ハンタウイルス vs その他ブニヤウイルス
- EEEV vs その他アルファウイルス（ゲタ、VEEV、WEEV）
- スプーマウイルス vs PERV（重要な識別）

### 5.3 データベース更新スケジュール

- **四半期レビュー:** NCBIで新規配列をチェック
- **年次再構築:** 最新リファレンスで全データベースを再構築
- **トリガー更新:** 臨床サンプルで新規変異株検出時

---- 

## 6. ストレージとバックアップ

### 6.1 プライマリストレージ
- **場所:** `/mnt/efs/databases/pmda/2024.1/`
- **ファイルシステム:** AWS EFS（伸縮性、EC2インスタンス間で共有）
- **アクセス:** 解析EC2インスタンス用読み取り専用
- **サイズ見積:**
  - ポリオーマウイルス：50 MB
  - ハンタウイルス：100 MB
  - アルファウイルス：200 MB
  - スプーマウイルス：50 MB
  - **合計:** \~400 MB（Kraken2インデックス含む：\~5 GB）

### 6.2 バックアップ戦略
- **S3バックアップ:** `s3://minion-data/databases/pmda/`
- **バックアップ頻度:** 週次自動同期
- **バージョニング:** 有効（S3バージョニング）
- **災害復旧:** S3から新規EFSマウントへ復元

---- 

## 7. データベース構築スクリプト

### 7.1 マスター構築スクリプト

```bash
#!/bin/bash
# build_pmda_4virus_databases.sh

set -e

BASE_DIR="/mnt/efs/databases/pmda/2024.1"
mkdir -p $BASE_DIR

# 1. ポリオーマウイルス
echo "ポリオーマウイルスデータベースを構築中..."
bash scripts/databases/build_polyomavirus_db.sh $BASE_DIR/polyomavirus

# 2. ハンタウイルス
echo "ハンタウイルスデータベースを構築中..."
bash scripts/databases/build_hantavirus_db.sh $BASE_DIR/hantavirus

# 3. アルファウイルス（EEEV）
echo "アルファウイルスデータベースを構築中..."
bash scripts/databases/build_alphavirus_db.sh $BASE_DIR/alphavirus

# 4. スプーマウイルス
echo "スプーマウイルスデータベースを構築中..."
bash scripts/databases/build_spumavirus_db.sh $BASE_DIR/spumavirus

# 5. S3へ同期
echo "S3へバックアップ中..."
aws s3 sync $BASE_DIR s3://minion-data/databases/pmda/2024.1/ --delete

echo "データベース構築完了！"
```

---- 

## 8. パイプラインとの統合

### 8.1 フェーズ4（病原体検出）統合

**現行フェーズ4データベース:**
- `/mnt/efs/databases/kraken2_standard/`
- `/mnt/efs/databases/blast_nt/`
- `/mnt/efs/databases/perv_references/`

**新規フェーズ4追加:**
```python
# scripts/phase4_pathogen/detect_pmda_4viruses.py

PMDA_4VIRUS_DATABASES = {
    "polyomavirus": "/mnt/efs/databases/pmda/2024.1/polyomavirus/polyoma_all.mmi",
    "hantavirus": "/mnt/efs/databases/pmda/2024.1/hantavirus/hantavirus_all.mmi",
    "alphavirus": "/mnt/efs/databases/pmda/2024.1/alphavirus/alphavirus_all.mmi",
    "spumavirus": "/mnt/efs/databases/pmda/2024.1/spumavirus/spumavirus_all_pol.fasta"  # BLASTのみ
}
```

---- 

## 9. 参考文献

1. NCBI Nucleotideデータベース: https://www.ncbi.nlm.nih.gov/nuccore/
2. NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy/
3. Minimap2: https://github.com/lh3/minimap2
4. Kraken2: https://github.com/DerrickWood/kraken2
5. BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi
6. プロトコル11: MinION\_Protocol\_11\_PMDA\_4ウイルス高感度検出プロトコル.md

