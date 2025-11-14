# 4-Virus Surveillance System

日本国内の畜産ブタにおける4種ウイルス（ハンタウイルス、ポリオーマウィルス、スピューマウィルス、東部ウマ脳炎ウイルス）の感染情報を統合監視するシステム

## システム概要

### 監視対象ウイルス
1. **ハンタウイルス (Hantavirus)** - 主にげっ歯類媒介
2. **ポリオーマウィルス (Polyomavirus)** - Sus scrofa polyomavirus
3. **スピューマウィルス (Spumavirus)** - MHLW特別管理病原体#5
4. **東部ウマ脳炎ウイルス (EEEV)** - 日本非流行

### 情報源
- **外部**:
  - MAFF（農林水産省サーベイランスレポート）
  - E-Stat（政府統計ポータル）
  - 学術論文（PubMed + J-STAGE）✅
- **内部**: MinIONパイプライン Phase 4検出結果（リアルタイム）

### アーキテクチャ
```
外部情報収集（日次）      内部パイプライン（リアルタイム）
     ↓                          ↓
  Lambda                    S3イベント
     ↓                          ↓
  DynamoDB  ←──────────────  Lambda
     ↓
 重要度分類エンジン
     ↓
 通知ルーター
     ├─ SNS/SES
     ├─ Slack ✅
     ├─ ダッシュボード
     └─ REST API
```

## ディレクトリ構造

```
surveillance/
├── external/              # 外部情報収集
│   ├── estat_client.py    # E-Stat APIクライアント
│   ├── maff_scraper.py    # MAFFスクレイパー
│   └── academic_monitor.py # 学術情報監視
├── internal/              # 内部パイプライン連携
│   ├── pipeline_listener.py # Phase 4結果監視
│   └── result_parser.py
├── alerting/              # 通知システム
│   ├── severity_engine.py      # 重要度分類
│   ├── notification_router.py  # 通知ルーティング
│   ├── slack_client.py         # Slack通知クライアント ✅
│   └── aws_notifier.py
├── dashboard/             # Streamlitダッシュボード
│   └── app.py
├── api/                   # REST API
│   └── main.py
├── lambda/                # Lambda関数
│   ├── external_collector/
│   ├── pipeline_listener/
│   └── alert_processor/
├── config/                # 設定ファイル
│   ├── severity_rules.yaml
│   └── config.yaml
└── models/                # データモデル
```

## セットアップ

### 1. 前提条件
- Python 3.11+
- AWS CLI設定済み（ap-northeast-1リージョン）
- Terraform 1.0+
- E-Stat APIアプリケーションID

### 2. 依存関係インストール
```bash
cd surveillance
pip install -r requirements.txt
```

### 3. 環境変数設定
```bash
# 基本設定
export SURVEILLANCE_BUCKET=surveillance-data
export E_STAT_APP_ID=bae1f981a6d093a9676b03c8eea37324b8de421b
export PUBMED_EMAIL=your-email@example.com
export AWS_REGION=ap-northeast-1

# Slack通知設定（オプション）
export SLACK_BOT_TOKEN=xoxb-your-token-here
export SLACK_APP_ID=A09TVLTGDSL
export SLACK_CLIENT_ID=580133244533.9947707557904
export SLACK_CLIENT_SECRET=your-client-secret
export SLACK_SIGNING_SECRET=your-signing-secret
```

または`.env`ファイルを使用:
```bash
cp surveillance/.env.template surveillance/.env
# .envファイルを編集してSlack認証情報を設定
export $(cat surveillance/.env | grep -v '^#' | xargs)
```

**Slack設定の詳細**: `surveillance/docs/SLACK_SETUP.md`を参照

### 4. AWSインフラデプロイ
```bash
cd infrastructure/surveillance
terraform init
terraform plan -var="estat_app_id=${E_STAT_APP_ID}" -var="pubmed_email=${PUBMED_EMAIL}"
terraform apply
```

### 5. DynamoDBテーブル作成（手動の場合）
```bash
aws dynamodb create-table --cli-input-json file://infrastructure/surveillance/dynamodb_schemas.json
```

## 使い方

### ダッシュボード起動
```bash
streamlit run surveillance/dashboard/app.py
```
→ http://localhost:8501

### REST API起動
```bash
uvicorn surveillance.api.main:app --reload --port 8000
```
→ http://localhost:8000/docs (API Documentation)

### 外部情報収集（手動実行）
```bash
# MAFF
python -m surveillance.external.maff_scraper

# E-Stat
python -m surveillance.external.estat_client

# Academic (PubMed + J-STAGE)
python surveillance/external/academic_monitor.py
# 出力例:
# === Testing PubMed Search ===
# PubMed: 5 articles found
# === Testing J-STAGE Search ===
# J-STAGE: 12 articles found
```

### Phase 4検出モジュール（既存パイプラインに統合）
```bash
./scripts/phase4_pathogen/detect_4viruses.py \
  --bam results/phase4/RUN-001/aligned.bam \
  --run-id RUN-001 \
  --sample-id SAMPLE-001 \
  --output-dir results/surveillance/
```

## API エンドポイント

### 検出情報
- `GET /api/v1/detections` - 検出一覧（フィルタ可能）
- `GET /api/v1/detections/{id}` - 特定検出の詳細
- `GET /api/v1/detections/realtime` - リアルタイム検出

### アラート
- `GET /api/v1/alerts/active` - アクティブアラート

### 外部情報源
- `GET /api/v1/external/daily-updates` - 日次更新情報

### 統計
- `GET /api/v1/statistics/trends` - トレンド分析
- `GET /api/v1/statistics/summary` - サマリー

## 重要度レベル

| レベル | 説明 | 対応時間 | 通知方法 |
|--------|------|----------|----------|
| **CRITICAL** | 即座対応必要（PERV同等） | < 5分 | SNS即時、SMS、Slack、ダッシュボード点滅 |
| **HIGH** | 緊急対応必要 | < 30分 | SNS、メール、Slack、ダッシュボード警告 |
| **MEDIUM** | 調査推奨 | < 2時間 | メール、Slack、ダッシュボード表示 |
| **LOW** | 情報記録 | < 24時間 | ダッシュボード記録 |

### Slackチャンネルルーティング
- **CRITICAL** → `#critical-alerts`
- **HIGH** → `#pathogen-alerts`
- **MEDIUM** → `#pathogen-monitoring`
- 日次サマリー → `#pathogen-monitoring`

### ウイルス別基準

#### スピューマウィルス
- `> 500 copies/mL` → CRITICAL
- `> 100 copies/mL` → HIGH
- `> 50 copies/mL` → MEDIUM
- 検出あり → LOW

#### ハンタウイルス
- `> 100 copies/mL` → HIGH（ブタでの検出は異常）
- 検出あり → LOW（げっ歯類汚染の可能性）

#### EEEV
- **任意の検出** → CRITICAL（日本非流行のため即座確認必要）

## 設定

### severity_rules.yaml
重要度分類ルールを定義：
- ウイルス別閾値
- 情報源別優先度
- 複合検出ルール
- 通知先設定

### config.yaml
システム全体設定：
- AWS設定
- API認証情報
- スケジュール設定
- データ保持期間

## 監視・運用

### CloudWatch Logs
```bash
aws logs tail /aws/lambda/surveillance-external-collector --follow
```

### DynamoDBクエリ例
```bash
# 最近のCRITICAL検出
aws dynamodb query \
  --table-name surveillance-detections \
  --index-name severity-index \
  --key-condition-expression "severity = :sev" \
  --expression-attribute-values '{":sev":{"S":"critical"}}'
```

### S3データ確認
```bash
aws s3 ls s3://surveillance-data/external/maff/2024/11/14/
```

## トラブルシューティング

### Lambda関数がタイムアウトする
- タイムアウト設定: 900秒（15分）まで延長可能
- メモリ: 512MB → 1024MBへ増量

### E-Stat APIエラー
- アプリケーションID確認: `config.yaml`のestat.app_id
- レート制限: 1秒あたり5リクエスト

### ダッシュボードが遅い
- キャッシュ設定: `config.yaml`のdashboard.cache_ttl調整
- データ範囲制限: 時間範囲を短縮

### 通知が届かない
- SNSトピックサブスクリプション確認
- SESメール検証ステータス確認
- Slack Bot Tokenが正しく設定されているか確認
- Slackチャンネルにbotが招待されているか確認
- CloudWatch Logsでエラー確認

## テスト

### 単体テスト
```bash
pytest tests/
```

### E-Stat API接続テスト
```bash
python -c "from surveillance.external.estat_client import EStatClient; client = EStatClient(); print(client.get_stats_list())"
```

### 重要度分類テスト
```bash
python surveillance/alerting/severity_engine.py
```

### Slack通知テスト
```bash
# 接続テスト
python surveillance/tests/test_slack_integration.py --test-conn

# アラート送信テスト
python surveillance/tests/test_slack_integration.py --test-alert

# 完全テストスイート
python surveillance/tests/test_slack_integration.py
```

## セキュリティ

- S3暗号化: AES256（デフォルト有効）
- DynamoDB暗号化: at-rest暗号化有効
- 機密情報: AWS Secrets Managerに保存推奨
- IAMロール: 最小権限原則

## パフォーマンス

- **外部収集**: 日次1回（11:00 JST）
- **内部監視**: リアルタイム（S3イベントトリガー）
- **ダッシュボード**: 30秒自動更新
- **API**: 最大1000件/リクエスト

## ライセンス

内部プロジェクトライセンス

## 連絡先

- **Critical Alert**: surveillance-critical@example.com
- **Technical Support**: support@example.com
- **Slack**: #surveillance-tech

## 変更履歴

### v2.2.0 (2025-11-15)
- ✅ **Slack通知統合実装**
  - Bot API + Webhook対応
  - 重要度別チャンネルルーティング
  - Block Kit形式のリッチメッセージ
  - アクションボタン（Critical アラート）
  - 日次サマリー通知
  - 包括的テストスイート
  - セットアップ自動化スクリプト
- 詳細ドキュメント: `surveillance/docs/SLACK_SETUP.md`

### v1.1.0 (2024-11-14)
- ✅ J-STAGE Webスクレイピング実装
- 日本語・英語検索対応
- PubMed + J-STAGE統合監視
- フォールバックパーシング機能

### v1.0.0 (2024-11-14)
- 初期リリース
- 外部情報収集（MAFF、E-Stat、PubMed）
- 内部パイプライン統合
- 重要度分類エンジン
- リアルタイムダッシュボード
- REST API
- AWS Lambda/DynamoDB/SNS統合

## 今後の拡張

- [x] J-STAGE Web スクレイピング統合 ✅
- [x] Slack通知実装 ✅
- [ ] 機械学習ベースのトレンド予測
- [ ] 地理的ヒートマップ
- [ ] モバイルアプリ対応
- [ ] 多言語対応（英語）
