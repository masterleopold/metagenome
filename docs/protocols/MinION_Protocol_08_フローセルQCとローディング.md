# MinION用メタゲノム解析プロトコル
# 第8章: フローセルQCとローディング

## 目次

1. [フローセルQCプロトコル概要](#1-フローセルqcプロトコル概要)
2. [MinION Platform QC実施手順](#2-minion-platform-qc実施手順)
3. [フローセルプライミング詳細プロトコル](#3-フローセルプライミング詳細プロトコル)
4. [ライブラリローディング計算と希釈](#4-ライブラリローディング計算と希釈)
5. [MinKNOWソフトウェア設定](#5-minknowソフトウェア設定)
6. [シーケンスラン開始とモニタリング](#6-シーケンスラン開始とモニタリング)
7. [リアルタイムQCメトリクス監視](#7-リアルタイムqcメトリクス監視)
8. [ラン時間最適化戦略](#8-ラン時間最適化戦略)
9. [トラブルシューティング](#9-トラブルシューティング)
10. [ALCOA+準拠記録](#10-alcoa準拠記録)

----

## 1. フローセルQCプロトコル概要

### 1.1 フローセルQCの目的

**MinION R9.4.1/R10.4.1フローセル品質確認の重要性**

フローセルQC（Platform QC）は、シーケンスラン開始前に実施する必須ステップです。このQCにより以下を確認します:

1. **アクティブポア数の確認**
   - 目標値: **1,200ポア以上** (新品フローセルの場合)
   - 最低許容値: **800ポア** (メタゲノム解析の場合)
   - 推奨値: **1,500ポア以上** (最適なスループット)

2. **フローセル状態の評価**
   - ポア劣化パターンの確認
   - 電流分布の異常検出
   - ナノポアメンブレンの完全性

3. **シーケンス予測スループットの推定**
   - 期待される総リード数
   - 目標カバレッジの達成可能性
   - ラン時間の最適化判断

**PMDA 91病原体検出におけるフローセルQC基準**

```yaml
フローセル品質基準:
  新品フローセル (推奨):
    アクティブポア数: ≥1,500
    予測収量: 15-30 Gb (48時間ラン)
    対応分析: DNA + RNA両ライブラリの同時解析

  使用可能フローセル (最低基準):
    アクティブポア数: ≥800
    予測収量: 8-15 Gb (48時間ラン)
    対応分析: DNA単独またはRNA単独

  使用不可フローセル:
    アクティブポア数: <800
    対応: 新しいフローセルへの交換必須
```

### 1.2 フローセルQC実施タイミング

**推奨実施スケジュール**

| タイミング | QC実施内容 | 判断基準 |
|----------|-----------|---------|
| フローセル受領時 | Platform QC | アクティブポア数確認、ロット記録 |
| ライブラリ調製完了後 | 再QC実施 | 保管期間が5日以上の場合 |
| シーケンス直前 (必須) | 最終QC | ローディング可否の最終判断 |
| ラン途中 (オプション) | ポア数モニタリング | ポア消費率の確認 |

**重要**: フローセルは使用まで4°Cで保管し、**受領後3ヶ月以内**に使用することを推奨します。

----

## 2. MinION Platform QC実施手順

### 2.1 必要機器・ソフトウェア

**ハードウェア要件**

```yaml
MinIONデバイス:
  モデル: MinION Mk1D または Mk1C
  USB接続: USB-C（Mk1Dの場合、USB 2.0転送速度対応）
  温度制御: ペルチェ素子搭載（Mk1D、10-35°C動作保証）
  電源: 独立電源アダプター使用推奨

コンピュータ仕様:
  OS: Ubuntu 20.04 LTS / Windows 10 Pro / macOS 10.15以降
  CPU: Intel Core i7以上 (8コア推奨)
  RAM: 16 GB以上 (32 GB推奨)
  ストレージ: 2 TB以上の空き容量 (NVMe SSD推奨)

環境条件:
  温度: 18-28°C (最適: 20-25°C)
  湿度: 30-70%
  振動: 最小限 (安定した実験台)
```

**ソフトウェア要件**

```yaml
MinKNOW:
  バージョン: 23.04以降 (最新版推奨)
  ベースコーラー: Dorado 0.5.0以降
  Duplexモード: 有効化必須

必須コンポーネント:
  - Guppy Basecaller (バックエンド)
  - MinKNOW Core Service
  - MinKNOW UI
  - Bream (Run Manager)
```

### 2.2 Platform QC実施プロトコル

**ステップバイステップ手順**

#### ステップ1: MinIONデバイスの準備

```bash
# 1. MinIONをUSB 3.0ポートに接続
#    - デバイスのLEDが青色に点灯することを確認
#    - MinKNOWが自動的にデバイスを認識

# 2. MinKNOWソフトウェアの起動
#    - Webブラウザで http://localhost:8000 にアクセス
#    - デバイス接続確認: "Device connected" 表示を確認

# 3. フローセルの取り出し
#    - フローセルを4°C冷蔵庫から取り出し
#    - 室温で30分間平衡化 (結露防止)
#    - パッケージを開封せずに平衡化
```

**重要な注意事項**:
- フローセルを急激に温度変化させない
- 開封前に外装の破損がないか確認
- シリアル番号を記録 (ALCOA+準拠)

#### ステップ2: フローセルのインストール

```yaml
フローセル装着手順:
  1. MinIONカバーをスライドして開く:
     - 右側のラッチを押して解除
     - カバーを右方向にスライド

  2. フローセルの取り出し:
     - 保護パッケージを開封
     - フローセルを静かに取り出す
     - フローセル底面のゴールドピンに触れない

  3. フローセルの装着:
     - MinIONの挿入スロットに合わせる
     - フローセルの向きを確認 (ポートが左側)
     - 静かに押し込み、カチッと音がするまで挿入

  4. プライミングポートキャップの確認:
     - スポットオンポートキャップが閉じていることを確認
     - プライミングポートキャップを時計回りに回して開く

  5. MinIONカバーの閉鎖:
     - カバーを左方向にスライド
     - ラッチがロックされることを確認
```

**チェックリスト**:
- [ ] フローセルが完全に挿入されている
- [ ] MinKNOWがフローセルを認識 (シリアル番号表示)
- [ ] プライミングポートキャップが開いている
- [ ] LEDが緑色に変化 (準備完了)

#### ステップ3: Platform QC実行

**MinKNOW UIでのQC開始**

```yaml
QC実行手順:
  1. MinKNOW UI画面で "Start" をクリック

  2. "Check flow cell" を選択:
     - Flow cell type: FLO-MIN106D (R9.4.1) または FLO-MIN114 (R10.4.1)
     - Kit selection: Platform QC

  3. Run nameの設定:
     - 命名規則: "PlatformQC_YYYYMMDD_FlowcellID"
     - 例: "PlatformQC_20250308_FAU12345"

  4. "Start platform QC" をクリック:
     - QC実行時間: 約10分
     - リアルタイムでポア数が表示される

  5. QC結果の確認:
     - Active pore count: 記録
     - Pore distribution: 可視化確認
     - Current histogram: 異常ピークの確認
```

**QC実行中のモニタリング**

MinKNOWは以下の情報をリアルタイム表示します:

```yaml
表示メトリクス:
  Active pores:
    - 総アクティブポア数
    - ポアグループ別分布 (1-512チャンネル)
    - 時間経過に伴うポア消費率

  Current distribution:
    - Open pore current (正常: 180-200 pA @ R9.4.1)
    - Blocked pore current (異常検出)
    - Adapter current (アダプター結合状態)

  Temperature:
    - フローセル温度 (目標: 34°C)
    - 温度安定性 (±0.5°C以内)
```

### 2.3 QC結果の判定基準

**アクティブポア数による判定**

```yaml
判定基準:
  優良 (Excellent):
    ポア数: ≥1,500
    対応: 高スループット解析に最適
    予測データ量: 20-30 Gb (48時間)
    推奨用途: DNA + RNA dual library解析

  良好 (Good):
    ポア数: 1,200-1,499
    対応: 標準的なメタゲノム解析に適用可能
    予測データ量: 15-20 Gb (48時間)
    推奨用途: DNA単独または RNA単独解析

  許容範囲 (Acceptable):
    ポア数: 800-1,199
    対応: 限定的な解析に使用可能
    予測データ量: 8-15 Gb (48時間)
    推奨用途: ターゲット病原体のみ (multiplexing制限)
    注意: RNA解析は非推奨 (収量不足リスク)

  不合格 (Failed):
    ポア数: <800
    対応: 使用不可、新しいフローセルへの交換
    理由: データ量不足により91病原体検出不可能
```

**ポア分布パターンによる判定**

正常なフローセルでは、512チャンネル全体に均一にポアが分布します:

```yaml
正常パターン:
  - 各チャンネルグループ (4グループ) で均等なポア数
  - ポア数の標準偏差: <15%
  - Dead zoneなし (全チャンネルで検出)

異常パターン (使用不可):
  - 特定チャンネルグループでポア数が50%以下
  - 連続した100チャンネル以上でポアゼロ
  - ポア数の標準偏差: >30%

要注意パターン (使用条件付き):
  - 1チャンネルグループのみポア数低下 (20-50%低下)
  - 対応: 低下グループを除外してラン実行
```

**電流ヒストグラムによる判定**

```yaml
Open pore current (R9.4.1フローセル):
  正常範囲: 180-200 pA
  ピーク位置: 190 pA付近
  ピーク幅: 狭く、シャープ

  異常パターン:
    - ピークが150 pA以下にシフト: ポア劣化
    - ピークが2つに分裂: 部分的ブロッキング
    - 幅広いピーク (>40 pA幅): 不安定なポア状態

Open pore current (R10.4.1フローセル):
  正常範囲: 140-160 pA
  ピーク位置: 150 pA付近
  ピーク幅: 狭く、シャープ
```

### 2.4 QC結果の記録

**ALCOA+準拠記録フォーマット**

```yaml
Platform QC記録シート:
  基本情報:
    - 実施日時: YYYY-MM-DD HH:MM (ISO 8601形式)
    - 実施者: フルネーム + 署名
    - フローセルID: FAU/FAT番号
    - ロット番号: フローセル箱に記載
    - 有効期限: YYYY-MM-DD

  QC結果:
    - Active pore count: 数値記録
    - Pore distribution: スクリーンショット添付
    - Current histogram: スクリーンショット添付
    - 判定結果: Excellent/Good/Acceptable/Failed

  環境条件:
    - 室温: °C
    - 湿度: %
    - MinKNOWバージョン: x.xx.x

  使用可否判定:
    - 使用可/不可: 判定理由記載
    - 推奨用途: DNA/RNA/両方
    - 予測データ量: Gb

  承認:
    - 確認者: フルネーム + 署名
    - 承認日時: YYYY-MM-DD HH:MM
```

**電子記録の保存**

```bash
# MinKNOW QCデータの保存場所
/var/lib/minknow/data/platform_qc/YYYYMMDD_HHMM_FAUxxxxx/

保存ファイル:
  - report_FAUxxxxx.md          # QCサマリーレポート
  - pore_activity.csv           # ポア数時系列データ
  - current_histogram.png       # 電流分布グラフ
  - channel_states.csv          # 全512チャンネル状態
  - duty_time.csv               # 各状態の時間割合
```

**記録の長期保存**: 全QCデータをバックアップサーバーに保存し、**最低5年間**保管します (PMDA規制対応)。

----

## 3. フローセルプライミング詳細プロトコル

### 3.1 プライミングの目的と原理

**プライミングの重要性**

フローセルプライミングは、シーケンスラン前に実施する必須ステップです:

1. **保存バッファーの除去**
   - フローセルは保存バッファー (Storage Buffer) で充填されて出荷
   - 保存バッファーはシーケンスに適さない組成
   - プライミングによりシーケンス用バッファーに置換

2. **ナノポアの活性化**
   - 適切なイオン濃度の確立
   - 電気伝導性の最適化
   - ポアの安定化

3. **気泡の除去**
   - フローセル内の微小気泡除去
   - 均一な流路確保
   - シーケンス品質の向上

**プライミングバッファー組成**

```yaml
Flush Buffer (FB):
  組成:
    - Flush Tether (FLT): フローセル膜安定化剤
    - Nuclease-free water: 希釈用
    - 塩類: イオン濃度調整

  調製比率:
    - FLT: 30 μL
    - Nuclease-free water: 970 μL
    - 総量: 1,000 μL (1 mL)

  保存条件:
    - 調製後すぐに使用 (推奨)
    - 短期保存: 4°C、24時間以内
    - 長期保存: 不可 (都度調製)
```

### 3.2 プライミングバッファー調製

**試薬の準備**

```yaml
必要試薬 (SQK-LSK114キットに含まれる):
  - FLT: Flush Tether (赤色キャップチューブ)
  - LFB: Long Fragment Buffer (緑色キャップチューブ)
  - SFB: Short Fragment Buffer (黄色キャップチューブ)
  - Nuclease-free water: 別途準備 (分子生物学グレード)

必要器具:
  - 1.5 mL DNase/RNase-freeチューブ: 2本
  - P1000ピペット + フィルターチップ
  - P200ピペット + フィルターチップ
  - ボルテックスミキサー
  - ミニ遠心機
```

**Flush Buffer調製プロトコル**

```yaml
ステップ1: FLTの準備
  - FLTチューブをボルテックス: 5秒間
  - 短時間スピンダウン: 3,000 rpm、5秒
  - 氷上で保管

ステップ2: Flush Bufferの調製
  1. 新しい1.5 mLチューブにラベル貼付:
     - "Flush Buffer - YYYYMMDD"

  2. Nuclease-free waterの添加:
     - P1000ピペットで970 μLを測り取る
     - 1.5 mLチューブに添加

  3. FLTの添加:
     - P200ピペットで30 μLを測り取る
     - 同じチューブに添加
     - 総量: 1,000 μL

  4. 混合:
     - ピペッティング10回 (泡立て厳禁)
     - または、チューブを反転5回
     - スピンダウン: 3,000 rpm、5秒

  5. 室温で平衡化:
     - 15分間放置
     - 使用直前まで遮光
```

**品質確認**

- 外観: 透明、均一 (沈殿物なし)
- pH: 8.0±0.2 (pH試験紙で確認可)
- 気泡: なし

### 3.3 プライミング実施プロトコル

**第1回プライミング (Initial Priming)**

```yaml
目的: 保存バッファーの完全除去

ステップ1: プライミングポートの準備
  1. MinIONカバーを開く
  2. プライミングポートキャップを時計回りに回して開く
  3. スポットオンポートキャップは閉じたまま

ステップ2: フローセルからのバッファー除去
  1. P1000ピペットに新しいフィルターチップを装着
  2. プライミングポートからゆっくり吸引:
     - 吸引量: 20-30 μL程度
     - 速度: ゆっくり (3秒以上かけて)
     - 注意: 気泡を吸い込まない

  3. 吸引したバッファーを廃棄:
     - 専用廃液容器に廃棄
     - ピペットチップも廃棄 (再使用厳禁)

ステップ3: Flush Bufferの注入
  1. 新しいP1000フィルターチップを装着
  2. Flush Bufferを800 μL吸引:
     - ピペット内に気泡が入らないよう注意
     - チップ先端を液面下に保つ

  3. プライミングポートから注入:
     - 注入速度: ゆっくり (30秒以上かけて)
     - 気泡が入らないよう細心の注意
     - 注入中に抵抗を感じたら一旦停止し、再開

  4. プライミングポートキャップを閉じる:
     - 反時計回りに回してしっかり閉める

ステップ4: 平衡化待機
  - MinIONカバーを閉じる
  - 5分間待機 (バッファーが全体に行き渡るまで)
  - MinKNOWでポア数をモニタリング (オプション)
```

**重要な注意事項**:

⚠️ **気泡混入の防止**
- 気泡が入るとポアが破壊される危険性
- ピペッティングは非常にゆっくり実施
- チップ先端は常に液面下に保つ
- 注入時は45度の角度で

⚠️ **過剰な圧力の回避**
- 注入時に強い抵抗を感じたら中止
- フローセル内圧が高すぎる可能性
- 数秒待ってから再開

⚠️ **クロスコンタミネーション防止**
- 全てのステップで新しいフィルターチップ使用
- ピペット本体をFlush Bufferに触れさせない

**第2回プライミング (Final Priming)**

ライブラリローディング直前に実施します (後述のセクション4で詳述)。

```yaml
目的: ライブラリローディング直前の最終準備

タイミング: ライブラリ調製完了後、ローディングの5分前

手順:
  1. プライミングポートから200 μL吸引
  2. Flush Bufferを200 μL注入
  3. プライミングポートキャップを開いたまま維持
     (ライブラリローディング時にそのまま使用)
```

### 3.4 プライミング後の確認

**MinKNOWでのポア数確認**

プライミング後、アクティブポア数が維持されているか確認します:

```yaml
確認項目:
  ポア数変動:
    - 許容範囲: ±5% (プライミング前後)
    - 例: 1,500ポア → 1,425-1,575ポア
    - 許容外: 新しいフローセルへ交換検討

  ポア分布:
    - 均一性維持の確認
    - 特定チャンネルグループの消失なし

  電流安定性:
    - Open pore currentピークの維持
    - ノイズレベルの増加なし
```

**プライミング不良の兆候と対処**

```yaml
症状1: ポア数が20%以上減少
  原因:
    - 気泡混入
    - 過剰な圧力
    - Flush Buffer不良
  対処:
    - 追加でFlush Buffer 200 μL注入
    - 10分間待機後、再度ポア数確認
    - 改善しない場合: フローセル交換

症状2: 電流ノイズの増加
  原因:
    - 不完全な混合
    - 温度不均一
  対処:
    - 5分間追加待機
    - MinION温度制御の確認 (34°C目標)

症状3: ポア分布の不均一化
  原因:
    - 注入速度が速すぎた
    - フローセル内圧の不均一
  対処:
    - ゆっくりとしたFlush Buffer追加注入 (100 μL)
    - 10分間平衡化
```

----

## 4. ライブラリローディング計算と希釈

### 4.1 ローディング量の決定

**最適ローディング量の計算**

MinIONフローセルへのライブラリローディング量は、アクティブポア数と目標データ量に基づいて決定します。

**基本計算式**

```yaml
ローディング量計算:
  基本式:
    Loading amount (fmol) = Active pores × Occupancy rate × Pore capacity

  パラメータ:
    Active pores: Platform QCで測定したポア数
    Occupancy rate: 30-50% (推奨: 40%)
    Pore capacity: 1ポアあたり1分子 (理論値)

  実用的な計算:
    Standard loading: 50-100 fmol (1,500ポア以上のフローセル)
    Low pore loading: 30-50 fmol (800-1,200ポアのフローセル)
```

**PMDA 91病原体解析における推奨ローディング量**

```yaml
DNA単独ライブラリ:
  新品フローセル (>1,500ポア):
    - ローディング量: 75 fmol
    - DNA量: 150-200 ng (平均サイズ10 kb想定)
    - 予測スループット: 20-25 Gb (48時間)

  使用済みフローセル (800-1,200ポア):
    - ローディング量: 40 fmol
    - DNA量: 80-120 ng
    - 予測スループット: 10-15 Gb (48時間)

RNA単独ライブラリ:
  新品フローセル (>1,500ポア):
    - ローディング量: 50 fmol
    - cDNA量: 100-150 ng (平均サイズ8 kb想定)
    - 予測スループット: 15-20 Gb (48時間)

  使用済みフローセル:
    - RNAライブラリは非推奨 (収量不足リスク)
```

### 4.2 DNA濃度からfmol変換

**モル濃度計算式**

```yaml
ステップ1: DNA濃度をng/μLからnM変換

  計算式:
    Molarity (nM) = (DNA concentration ng/μL × 10^6) / (Average size bp × 650)

  定数:
    - 10^6: ng → g 変換 + nM単位変換
    - 650: 1 bp あたりの平均分子量 (g/mol)

  例: 20 ng/μL、平均サイズ10 kb (10,000 bp)
    Molarity = (20 × 10^6) / (10,000 × 650)
             = 20,000,000 / 6,500,000
             = 3.08 nM
             = 3,080 pM
             = 3.08 fmol/μL

ステップ2: 必要体積の計算

  計算式:
    Volume (μL) = Target fmol / Molarity (fmol/μL)

  例: 75 fmol ローディング、3.08 fmol/μL ライブラリ
    Volume = 75 / 3.08 = 24.4 μL
```

**実用的な計算表**

| DNA濃度 (ng/μL) | 平均サイズ (kb) | モル濃度 (fmol/μL) | 75 fmol必要量 (μL) | 50 fmol必要量 (μL) |
|----------------|----------------|-------------------|-------------------|-------------------|
| 15 | 8 | 2.88 | 26.0 | 17.4 |
| 20 | 8 | 3.85 | 19.5 | 13.0 |
| 25 | 8 | 4.81 | 15.6 | 10.4 |
| 15 | 10 | 2.31 | 32.5 | 21.7 |
| 20 | 10 | 3.08 | 24.4 | 16.3 |
| 25 | 10 | 3.85 | 19.5 | 13.0 |
| 15 | 12 | 1.92 | 39.0 | 26.0 |
| 20 | 12 | 2.56 | 29.3 | 19.5 |
| 25 | 12 | 3.21 | 23.4 | 15.6 |

### 4.3 ライブラリ希釈プロトコル

**希釈の必要性判定**

```yaml
希釈不要ケース:
  - ライブラリ濃度: <5 fmol/μL
  - 75 fmol必要量が15 μL以上
  - そのままローディング可能

希釈必要ケース:
  - ライブラリ濃度: >5 fmol/μL
  - 75 fmol必要量が10 μL未満
  - 理由: 小容量ではピペッティング誤差大
  - 対応: 希釈して15-30 μL使用
```

**希釈プロトコル**

```yaml
使用バッファー: Elution Buffer (EB) またはNuclease-free water

希釈計算:
  目標: 最終濃度 3-5 fmol/μL

  例: 現在15 ng/μL (10 kb平均) = 2.31 fmol/μL
    → 希釈不要 (そのまま32.5 μL使用)

  例: 現在30 ng/μL (10 kb平均) = 4.62 fmol/μL
    → 希釈不要 (そのまま16.2 μL使用)

  例: 現在50 ng/μL (10 kb平均) = 7.69 fmol/μL
    → 希釈必要
    希釈比率: 7.69 / 4.0 = 1.92倍
    調製例: ライブラリ10 μL + EB 9.2 μL = 4.0 fmol/μL
    必要量: 75 / 4.0 = 18.75 μL

希釈手順:
  1. 新しい0.2 mL PCRチューブを準備
  2. ラベル貼付: "Diluted Library - YYYYMMDD"
  3. EB/水を先に添加 (計算量)
  4. ライブラリを添加 (計算量)
  5. ピペッティング混合: 10回
  6. スピンダウン: 3,000 rpm、5秒
  7. 氷上で保管 (使用まで)
```

### 4.4 ライブラリミックス調製

**最終ローディングミックス組成**

MinIONへのローディングは、ライブラリを以下のバッファーと混合して実施します:

```yaml
SQK-LSK114 Loading Mix:

  構成要素:
    - Sequencing Buffer II (SB II): 37.5 μL
    - Loading Beads II (LB II): 25.5 μL
    - ライブラリ (DNA): 12 μL (例: 75 fmol)
    - 総量: 75 μL

  各成分の役割:
    Sequencing Buffer II:
      - イオン濃度調整
      - pH緩衝 (pH 8.0)
      - ナノポア安定化

    Loading Beads II:
      - ライブラリ分散
      - フローセル内均一分布促進
      - ポアへの効率的な供給
      - 重要: 使用直前にボルテックス必須

    ライブラリ:
      - アダプター付加済みDNA/cDNA
      - 濃度調整済み (3-5 fmol/μL)
```

**ローディングミックス調製プロトコル**

```yaml
ステップ1: 試薬の準備
  1. Sequencing Buffer II (SB II)を冷蔵庫から取り出し:
     - 室温で10分間平衡化
     - ボルテックス: 5秒間
     - スピンダウン: 3,000 rpm、5秒

  2. Loading Beads II (LB II)を準備:
     - ボルテックス: 30秒間 (強く振る)
     - 重要: ビーズが完全に懸濁するまで
     - スピンダウン: 3,000 rpm、5秒

  3. ライブラリを氷上から取り出し:
     - 室温で2分間平衡化
     - ピペッティング混合: 5回
     - スピンダウン: 3,000 rpm、5秒

ステップ2: ローディングミックスの調製
  1. 新しい1.5 mL DNA LoBind tubeを準備
  2. ラベル貼付: "Loading Mix - YYYYMMDD - Sample ID"

  3. 添加順序 (必ず守る):
     a. Sequencing Buffer II: 37.5 μL
        - P200ピペットで正確に測定

     b. ライブラリ: 12 μL
        - P20ピペットで正確に測定
        - チップ先端をバッファーに触れさせて排出

     c. Loading Beads II: 25.5 μL (最後に添加)
        - P200ピペットで正確に測定
        - 添加直前に再度チューブを反転混合

  4. 混合:
     - ピペッティング混合: 10回
     - 速度: ゆっくり (泡立て厳禁)
     - スピンダウン: 3,000 rpm、5秒

  5. 氷上で保管:
     - 5分以内にローディング実施
     - 長時間放置厳禁 (ビーズ沈降)
```

**チェックリスト**

ローディング直前確認事項:
- [ ] Loading Beads IIが完全に懸濁している
- [ ] 総量が75 μLである
- [ ] 気泡が混入していない
- [ ] ローディングミックス調製から5分以内
- [ ] プライミングポートが開いている
- [ ] MinKNOWソフトウェアがシーケンス準備完了

----

## 5. MinKNOWソフトウェア設定

### 5.1 MinKNOWシーケンス実験設定

**シーケンス実験の開始**

```yaml
MinKNOW UI操作手順:

ステップ1: "Start" ボタンをクリック
  - MinKNOW UIトップ画面
  - フローセルが接続されていることを確認

ステップ2: "Start sequencing" を選択
  - Platform QCではなくSequencingを選択

ステップ3: Kit selectionの設定
  Sequencing kit:
    - "SQK-LSK114" を選択
    - またはKit検索窓で "LSK114" と入力

  Expansion kit (オプション):
    - バーコード使用の場合: "EXP-NBD196" など選択
    - バーコード不使用: "No expansion kit" 選択

ステップ4: Experiment settingsの設定
  Run name:
    - 命名規則: "ProjectName_SampleID_YYYYMMDD_RunNumber"
    - 例: "Metagenome_PigPlasma001_20250308_Run01"
    - 重要: アンダースコア使用、スペース不可

  Sample ID:
    - サンプル固有ID入力
    - 例: "DonorPig_12345_Plasma_DNA"
    - トレーサビリティ確保

  Experiment group (オプション):
    - プロジェクトグループ名
    - 例: "PMDA_91Pathogen_Screening"
```

### 5.2 Basecalling設定 (Duplex Mode)

**Duplexベースコーリングの重要性**

PMDA 91病原体検出では、**Duplex basecalling**が必須です:

```yaml
Duplexモードの利点:
  精度向上:
    - Simplex: Q20-Q25 (99-99.9% accuracy)
    - Duplex: Q30以上 (99.9% accuracy)
    - 重要: 偽陽性病原体検出の劇的減少

  病原体同定の信頼性:
    - SNP level変異検出: 可能
    - 系統解析: 高精度
    - 耐性�り子検出: 信頼性向上

  データ効率:
    - 必要リード数: 1/10に減少
    - コスト削減: シーケンス深度要求低下
```

**Basecalling設定詳細**

```yaml
MinKNOW Basecalling Configuration:

ステップ1: Basecalling mode選択
  位置: "Output" タブ内

  選択肢:
    ○ Super accuracy (Duplex) - 推奨
      - Duplexベースコール有効
      - 精度: Q30以上
      - 速度: リアルタイム処理可能 (GPUあり)

    ○ High accuracy
      - Simplex basecallのみ
      - 精度: Q20-25
      - 使用ケース: GPU不足時の一時対応

    ○ Fast basecalling
      - 低精度
      - 病原体解析には不適
      - 使用非推奨

ステップ2: Basecaller configuration
  Model selection:
    - "dna_r10.4.1_e8.2_400bps_sup.cfg" (R10.4.1フローセル)
    - "dna_r9.4.1_450bps_sup.cfg" (R9.4.1フローセル)
    - 自動選択: MinKNOWが通常自動で正しいモデル選択

  Duplex settings:
    - Enable duplex basecalling: ☑ チェック
    - Duplex mode: "Automatic"
    - Simplex filtering: "Duplex pairs only" (オプション)
      → Duplex readsのみ出力 (データ量削減)

ステップ3: Output options
  Output folder:
    - デフォルト: /var/lib/minknow/data/
    - 推奨: 外部NASまたは大容量SSD
    - 必要容量: 1 TB以上 (48時間ラン)

  File format:
    - ☑ FASTQ (必須)
    - ☑ BAM (オプション: アライメント解析用)
    - ☐ FAST5 (オプション: 生シグナル保存、容量大)

  Read splitting:
    - "4000" reads per file (推奨)
    - 理由: ファイルサイズ管理、並列処理容易化
```

### 5.3 GPU設定とリソース管理

**GPU basecallingの重要性**

```yaml
GPU要件:
  推奨GPU:
    - NVIDIA RTX 4090 (24 GB VRAM)
    - NVIDIA A6000 (48 GB VRAM)
    - NVIDIA A100 (40/80 GB VRAM)

  最小GPU:
    - NVIDIA RTX 3070 (8 GB VRAM)
    - 性能: リアルタイム処理ギリギリ
    - Duplex処理遅延の可能性

  GPU不足時の対処:
    - オフラインbasecalling (ラン後処理)
    - Dorado standalone使用
    - クラウドGPU (AWS/Google Cloud)
```

**MinKNOW GPU設定**

```yaml
GPU configuration (MinKNOW):

  設定ファイル編集:
    - 場所: /opt/ont/minknow/conf/app_conf
    - ファイル: guppy_config.toml

  重要パラメータ:
    gpu_runners_per_device = 4
      - GPU並列処理数
      - 推奨: 4-8 (GPU VRAMに依存)

    chunks_per_runner = 512
      - 各ランナーの処理単位
      - 推奨: 512-1024

    chunk_size = 2000
      - シグナルチャンクサイズ
      - 推奨: 2000 (デフォルト)

  適用方法:
    sudo systemctl restart minknow
```

### 5.4 バーコード設定 (Multiplexing)

**24サンプル同時解析のためのバーコーディング**

```yaml
Native Barcode Kit使用:
  キット: EXP-NBD196 (96 barcodes)

  PMDA病原体スクリーニングでの使用例:
    - 1フローセル: 24サンプル multiplexing
    - 各サンプル期待収量: 0.8-1.2 Gb (48時間ラン)
    - 十分な深度: 病原体検出に対応

MinKNOW barcode設定:
  ステップ1: Expansion kit選択
    - "EXP-NBD196" を選択

  ステップ2: Barcode配置指定
    - "Specify barcodes" を選択
    - 使用するバーコード番号を入力
    - 例: BC01, BC02, ..., BC24

  ステップ3: Sample ID割り当て
    - 各バーコードに対応するサンプルIDを入力
    - CSV形式でのアップロードも可能

  ステップ4: Demultiplexing設定
    - Real-time demultiplexing: ☑ 有効化
    - Min score: 60 (推奨)
    - Require both ends: ☑ チェック (精度向上)

  出力:
    - バーコードごとに別フォルダ生成
    - 例: /fastq_pass/barcode01/, barcode02/, ...
```

### 5.5 Run duration設定

**ラン時間の最適化**

```yaml
Run duration options:

  推奨設定: Time-based run
    Duration: 48 hours
    理由:
      - 91病原体検出に十分なデータ量
      - Duplex read生成に必要な時間
      - ポア消費とデータ量のバランス最適

  代替設定: Target yield
    Target: 20 Gb
    理由:
      - データ量重視
      - 早期完了可能性
      - フローセル再利用可能

  緊急時設定: 6-12 hours
    - 緊急病原体スクリーニング
    - ターゲット病原体限定検出
    - データ量: 5-10 Gb

MinKNOW設定:
  位置: "Run options" タブ

  Time-based run:
    - Run duration: 48 hours
    - Auto-stop: ☑ 有効化

  または Yield-based run:
    - Target yield: 20 Gb
    - Auto-stop when reached: ☑ 有効化
```

----

## 6. シーケンスラン開始とモニタリング

### 6.1 ライブラリローディング実施

**ローディング直前最終確認**

```yaml
チェックリスト:
  MinKNOW設定:
    - [ ] Sequencing kit: SQK-LSK114
    - [ ] Basecalling: Duplex mode enabled
    - [ ] Run name: 正しく入力
    - [ ] Run duration: 48 hours設定

  ライブラリ準備:
    - [ ] Loading mix調製完了 (75 μL)
    - [ ] 調製から5分以内
    - [ ] 気泡なし
    - [ ] 氷上保管中

  フローセル状態:
    - [ ] プライミング完了
    - [ ] アクティブポア数: >1,200
    - [ ] プライミングポートキャップ: 開いている
    - [ ] スポットオンポートキャップ: 閉じている
```

**ローディング実施手順**

```yaml
ステップ1: Loading mixの最終準備
  1. Loading mixチューブをピペッティング混合:
     - 5回、ゆっくりと (ビーズ再懸濁)
     - 泡立て厳禁

  2. スピンダウン:
     - 3,000 rpm、5秒
     - 全液体をチューブ底に集める

ステップ2: プライミングポートからのバッファー除去
  1. MinIONカバーを開く
  2. プライミングポートキャップが開いていることを確認

  3. P1000ピペット + 新しいフィルターチップ装着
  4. プライミングポートから200 μL吸引:
     - ゆっくり (5秒以上かけて)
     - 気泡を吸い込まない

  5. 吸引した液体を廃棄

ステップ3: Loading mixの注入
  1. 新しいP1000フィルターチップを装着
  2. Loading mix全量 (75 μL) を吸引:
     - チューブ底からゆっくり吸引
     - ビーズが均一に吸引されるよう注意
     - 気泡混入厳禁

  3. プライミングポートから注入:
     - 注入速度: 非常にゆっくり (30秒以上)
     - ピペットを45度の角度で保持
     - 抵抗を感じたら一時停止し再開
     - 全量注入を確認

  4. プライミングポートキャップを閉じる:
     - 反時計回りに回してしっかり閉める

ステップ4: SpotONポートからの追加ローディング (オプション)
  ※ 通常はプライミングポートのみで十分
  ※ 高スループット要求時のみ実施

  1. スポットオンポートキャップを開く
  2. 残りのLoading mix (もし調製していれば) を滴下:
     - 1滴ずつ、ゆっくりと
     - 総量: 最大200 μL追加可能

  3. スポットオンポートキャップを閉じる

ステップ5: MinIONカバーを閉じる
  1. カバーを左方向にスライド
  2. ラッチがロックされることを確認
  3. LED表示確認: 緑色 (シーケンス準備完了)
```

**ローディング後の待機時間**

```yaml
平衡化時間:
  最小: 5分
  推奨: 10分

  目的:
    - ライブラリのフローセル全体への拡散
    - ポアへのDNA分子の結合
    - 電流の安定化

  この間の作業:
    - MinKNOWでシーケンス設定の最終確認
    - ローディング記録の記載 (ALCOA+準拠)
    - 環境条件の記録 (温度、湿度)
```

### 6.2 シーケンスラン開始

**MinKNOWでのラン開始**

```yaml
ステップ1: 設定の最終確認
  MinKNOW UI "Start" 画面で確認:
    - Kit: SQK-LSK114
    - Basecalling: Super accuracy (Duplex)
    - Run name: 正しいか
    - Run duration: 48 hours
    - Output folder: 容量十分か

ステップ2: "Start run" ボタンをクリック
  - クリック後、即座にシーケンス開始
  - 初期化: 約2分間
  - Basecalling開始: 初期化完了後

ステップ3: Run start確認
  確認項目:
    - Run status: "Running"
    - Active channels: 増加中
    - Reads/sec: >50 (安定後)
    - Basecalling queue: 増加中
```

**ラン開始直後のモニタリング (最初の30分)**

```yaml
重要メトリクス:

  1. Active channels (時間経過):
     0-5分: 100-300チャンネル
     5-15分: 500-1,000チャンネル
     15-30分: 1,000-1,500チャンネル (安定)

  2. Reads per second:
     0-5分: 10-30 reads/sec
     5-15分: 30-100 reads/sec
     15-30分: 100-200 reads/sec (安定)

  3. Basecalling speed:
     - Samples/sec: >2,000,000
     - GPU utilization: 70-90%
     - Queue length: <1,000 reads

  4. Read N50:
     - 初期: 5-10 kb
     - 30分後: 8-15 kb
     - 安定期: 10-20 kb

正常なラン開始の兆候:
  - [ ] 15分以内にActive channels >500
  - [ ] Reads/secが順調に増加
  - [ ] エラーメッセージなし
  - [ ] GPU basecallingがリアルタイム処理
  - [ ] Read N50 >5 kb

異常の兆候と対処:
  Active channels <200 (15分後):
    原因: ライブラリ濃度不足、ローディング失敗
    対処: 追加ライブラリローディング (50 fmol)

  Reads/sec <20 (15分後):
    原因: フローセル不良、プライミング不良
    対処: ラン中止、フローセル交換

  Basecalling queue >5,000:
    原因: GPU性能不足
    対処: オフラインbasecallingに切り替え
```

### 6.3 リアルタイムデータモニタリング

**MinKNOW UIダッシュボード**

MinKNOWは以下のリアルタイム情報を提供します:

```yaml
メイン画面:
  1. Run summary:
     - Elapsed time: 経過時間
     - Estimated end time: 終了予定時刻
     - Active channels: 現在のアクティブチャンネル数
     - Total reads: 総リード数
     - Total yield: 総データ量 (Gb)

  2. Channel states (512チャンネル可視化):
     - 緑: Sequencing (シーケンス中)
     - 青: Available (利用可能、DNA待機中)
     - 灰色: Inactive (不活性)
     - 赤: Saturated (過負荷)

  3. Read length histogram:
     - リード長分布
     - N50値
     - 最大リード長

  4. Duty time plot:
     - Sequencing時間割合
     - Adapter時間割合
     - Pore時間割合
```

**詳細メトリクス (Detailsタブ)**

```yaml
Throughput metrics:
  - Reads per second: リアルタイムスループット
  - Bases per second: 塩基配列生成速度
  - Instantaneous rate: 瞬間速度
  - Average rate: 平均速度

Quality metrics:
  - Mean Q score: 平均品質スコア
  - Median Q score: 中央値品質スコア
  - Q20+ pass rate: Q20以上の割合
  - Q30+ pass rate: Q30以上の割合 (Duplex)

Duplex metrics (Duplex mode時):
  - Duplex reads: Duplexリード数
  - Duplex rate: Duplex化率 (目標: >40%)
  - Simplex reads: Simplexリード数
  - Duplex yield: Duplexデータ量
```

----

## 7. リアルタイムQCメトリクス監視

### 7.1 時間経過ごとの正常値範囲

**ラン開始後1時間**

```yaml
期待されるメトリクス:

  Active channels:
    正常範囲: 1,000-1,500 (新品フローセル)
    最低許容: 800

  Total reads:
    正常範囲: 50,000-150,000 reads
    Read N50: 8-12 kb

  Total yield:
    正常範囲: 0.4-1.0 Gb
    平均Q score: Q15-20 (Simplexが主体)

  Duty time:
    Sequencing: 60-75%
    Adapter: 10-15%
    Pore: 10-20%

  Basecalling:
    Queue: <500 reads (GPU処理が追いついている)
    Samples/sec: >2,000,000
```

**ラン開始後6時間**

```yaml
期待されるメトリクス:

  Active channels:
    正常範囲: 900-1,400
    ポア消費率: 5-10%/6時間

  Total reads:
    正常範囲: 300,000-800,000 reads
    Read N50: 10-15 kb

  Total yield:
    正常範囲: 3-6 Gb
    Q30+ reads (Duplex): 開始増加
    Duplex rate: 20-40%

  Duty time:
    Sequencing: 65-80%
    Adapter: 8-12%
    Pore: 10-15%
```

**ラン開始後24時間**

```yaml
期待されるメトリクス:

  Active channels:
    正常範囲: 700-1,200
    ポア消費率: 20-30%

  Total reads:
    正常範囲: 1,000,000-2,500,000 reads
    Read N50: 12-18 kb

  Total yield:
    正常範囲: 10-18 Gb
    Q30+ reads: 増加中
    Duplex rate: 40-60%

  Duty time:
    Sequencing: 60-75%
    Adapter: 8-12%
    Pore: 10-20%
```

**ラン終了時 (48時間)**

```yaml
最終メトリクス目標:

  Active channels:
    正常範囲: 500-1,000
    総ポア消費: 40-60%

  Total reads:
    目標: 2,000,000-5,000,000 reads
    Read N50: 12-20 kb

  Total yield:
    目標: 20-30 Gb
    うちDuplex: 8-15 Gb (40-50%)

  Quality:
    Mean Q score: Q25-30
    Q30+ rate: 40-60% (Duplex込み)

  PMDA病原体検出に対する評価:
    合格基準:
      - Total yield: ≥15 Gb
      - Q30+ yield: ≥5 Gb
      - Read N50: ≥10 kb
      - Active channels (終了時): ≥500
```

### 7.2 アラート設定と対応

**MinKNOW自動アラート**

```yaml
Critical alerts (即座対応必要):

  1. "Low active channels" (アクティブチャンネル急減)
     条件: Active channels <200
     原因:
       - フローセル不良
       - 温度制御失敗
       - 電源トラブル
     対応:
       - ラン中止
       - エラーログ確認
       - フローセル交換

  2. "Basecalling queue overflow"
     条件: Queue >10,000 reads
     原因: GPU性能不足
     対応:
       - オフラインbasecallingに切り替え
       - GPU設定見直し

  3. "Temperature control error"
     条件: フローセル温度 <30°C または >38°C
     原因:
       - MinION冷却ファン故障
       - 室温異常
     対応:
       - 室温確認、エアコン調整
       - MinION再起動
       - 改善しない場合: ラン中止

Warning alerts (監視強化):

  1. "Low sequencing duty time"
     条件: Sequencing duty <50%
     原因:
       - ライブラリ濃度不足
       - アダプター不良
     対応:
       - 追加ライブラリローディング検討
       - 次回ライブラリ調製の改善

  2. "Low N50"
     条件: Read N50 <5 kb
     原因:
       - DNA断片化
       - サイズ選択不良
     対応:
       - データ取得継続 (短鎖リードでも情報あり)
       - 次回プロトコル改善

  3. "Low duplex rate"
     条件: Duplex rate <20% (24時間後)
     原因:
       - シーケンス深度不足
       - DNA鎖対形成不足
     対応:
       - ラン時間延長検討
       - 次回ライブラリ濃度最適化
```

### 7.3 中間データ確認

**FASTQ出力ファイルの確認**

ラン進行中でも、FASTQ fileは逐次生成されます:

```yaml
出力ディレクトリ:
  /var/lib/minknow/data/[run_name]/fastq_pass/

ファイル命名規則:
  FAUxxxxx_pass_barcode01_[file_number].fastq.gz

中間確認タイミング:
  - 6時間後
  - 24時間後

確認コマンド例:

  # FASTQ file数の確認
  ls -lh /var/lib/minknow/data/[run_name]/fastq_pass/ | wc -l

  # 総リード数の概算
  zcat /var/lib/minknow/data/[run_name]/fastq_pass/*.fastq.gz | \
    grep -c "^@"

  # 平均リード長の確認 (最初の10,000リード)
  zcat /var/lib/minknow/data/[run_name]/fastq_pass/*.fastq.gz | \
    head -40000 | \
    awk 'NR%4==2{sum+=length($0); count++} END{print sum/count}'

  # Q30リードの割合確認
  zcat /var/lib/minknow/data/[run_name]/fastq_pass/*.fastq.gz | \
    NanoStat --fastq stdin | grep "Mean read quality"
```

**簡易病原体検出テスト (オプション)**

24時間時点で、簡易的な病原体検出テストを実施できます:

```bash
# Kraken2による高速分類 (中間FASTQ使用)
kraken2 --db /path/to/kraken2_db \
        --threads 16 \
        --output kraken2_intermediate.txt \
        --report kraken2_report_intermediate.txt \
        /var/lib/minknow/data/[run_name]/fastq_pass/*.fastq.gz

# レポート確認
head -50 kraken2_report_intermediate.txt

# 主要病原体の検出確認
grep "Porcine circovirus" kraken2_report_intermediate.txt
grep "Pseudorabies virus" kraken2_report_intermediate.txt
```

このテストにより、ラン継続の妥当性を判断できます。

----

## 8. ラン時間最適化戦略

### 8.1 目標データ量による早期終了判断

**PMDA病原体検出に必要な最小データ量**

```yaml
DNA単独ライブラリ:
  最小要求: 10 Gb (Q20以上)
  推奨: 15-20 Gb
  最適: 25-30 Gb

  早期終了判断:
    20 Gb達成 + Read N50 >10 kb + Q30 rate >30%
    → 24-30時間で達成可能な場合、早期終了検討

RNA単独ライブラリ (cDNA):
  最小要求: 8 Gb (Q20以上)
  推奨: 12-18 Gb
  最適: 20-25 Gb

  早期終了判断:
    15 Gb達成 + Read N50 >8 kb + Q30 rate >25%
    → 20-30時間で達成可能

Multiplex (24サンプル):
  最小要求: 24 Gb (各サンプル1 Gb)
  推奨: 30-36 Gb (各サンプル1.5 Gb)

  早期終了判断:
    全サンプルで1.5 Gb以上 + バランス良好
    → 30-40時間で達成可能
```

**早期終了の実施手順**

```yaml
MinKNOWでの停止:
  1. MinKNOW UI画面で "Stop run" をクリック
  2. 理由選択: "Target yield reached"
  3. 確認ダイアログ: "Stop run" をクリック

  停止処理:
    - シーケンス即座停止
    - Basecalling queue処理完了まで継続
    - 最終レポート生成
    - 所要時間: 5-15分

  データ保存:
    - 全FASTQファイル保存済み
    - Summary report生成済み
    - MinKNOWデータベースに記録
```

### 8.2 フローセル再利用の判断

**再利用可能性の評価**

```yaml
再利用推奨条件:
  - ラン終了時アクティブポア数: ≥800
  - ラン時間: ≤24時間
  - 総データ量: ≤15 Gb
  - ポア消費率: ≤40%

再利用手順:
  1. ラン終了後、Nuclease Flush実施:
     - Flow Cell Wash Kit (EXP-WSH004) 使用
     - プロトコル時間: 30分

  2. 再度Platform QC実施:
     - アクティブポア数確認
     - 800ポア以上なら再利用可能

  3. 新しいライブラリでシーケンス再開:
     - プライミングから再実施
     - 期待収量: 初回の60-80%

再利用非推奨条件:
  - アクティブポア数: <600
  - 電流分布異常
  - 48時間フルラン実施済み
  → 新しいフローセルへの交換推奨
```

### 8.3 ラン時間とコストのバランス

**ラン時間別のコストパフォーマンス**

| ラン時間 | 期待データ量 | フローセルあたりコスト | Gb単価 | 推奨用途 |
|---------|------------|---------------------|-------|---------|
| 6時間 | 3-5 Gb | ¥80,000 | ¥16,000-26,000 | 緊急スクリーニング |
| 12時間 | 6-10 Gb | ¥80,000 | ¥8,000-13,000 | 迅速検査 |
| 24時間 | 12-18 Gb | ¥80,000 | ¥4,400-6,700 | 標準検査 |
| 48時間 | 20-30 Gb | ¥80,000 | ¥2,700-4,000 | 最適 (推奨) |
| 72時間 | 25-35 Gb | ¥80,000 | ¥2,300-3,200 | 延長 (ポア消費大) |

**最適化戦略**

```yaml
コスト重視:
  - ラン時間: 48-72時間
  - データ量最大化
  - Gb単価最小化
  - 推奨: 通常スクリーニング

速度重視:
  - ラン時間: 6-12時間
  - 緊急病原体検出
  - 高コストだが迅速
  - 推奨: アウトブレイク対応

バランス型:
  - ラン時間: 24-36時間
  - 十分なデータ量
  - 合理的コスト
  - 推奨: 定期スクリーニング
```

----

## 9. トラブルシューティング

### 9.1 低スループット問題

**症状: Total yield <5 Gb (24時間後)**

```yaml
原因1: ライブラリ濃度不足
  確認方法:
    - Active channels <500
    - Sequencing duty time <40%

  対処法:
    a. 追加ライブラリローディング:
       - 新たに50 fmol準備
       - Loading mix調製 (前述プロトコル)
       - プライミングポートから注入

    b. 次回改善:
       - ライブラリ調製でのDNA回収率向上
       - サイズ選択最適化
       - Qubit濃度測定の精度確認

原因2: フローセル不良
  確認方法:
    - 初回Platform QCでポア数が低かった (<1,000)
    - ポア消費が異常に速い (>50%/24時間)

  対処法:
    - ラン中止
    - 新しいフローセルへ交換
    - 不良フローセルを記録、Oxford Nanoporeへ報告

原因3: プライミング不良
  確認方法:
    - Channel states画面で不均一なポア分布
    - 特定領域でポアゼロ

  対処法:
    - 追加プライミング実施:
      1. ラン一時停止
      2. Flush Buffer 500 μL注入
      3. 10分間待機
      4. ラン再開
```

### 9.2 低Duplex率問題

**症状: Duplex rate <20% (24時間後)**

```yaml
原因1: シーケンス深度不足
  確認方法:
    - Total reads <500,000 (24時間後)
    - Reads/sec <50

  対処法:
    - ラン時間延長 (48時間→72時間)
    - 追加ライブラリローディング

原因2: DNAフラグメントサイズ不適
  確認方法:
    - Read N50 <5 kb
    - リード長分布が短鎖に偏る

  対処法 (次回改善):
    - AMPure XP比率変更: 0.4× → 0.6×
    - DNA抽出時のシェアリング最小化
    - ピペッティング速度をさらに緩やかに

原因3: Basecalling設定ミス
  確認方法:
    - MinKNOW設定確認: Duplex mode OFF?

  対処法:
    - ラン停止
    - オフラインDuplex basecalling実施:

      # Dorado duplex basecalling
      dorado duplex \
        dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
        /path/to/pod5_files/ \
        --emit-fastq > duplex_reads.fastq
```

### 9.3 ポア急速消費問題

**症状: Active channels <500 (12時間後)**

```yaml
原因1: ライブラリ過剰ローディング
  確認方法:
    - 初期Active channels >1,500
    - 6時間後に急減 (<800)
    - Channel states: 多数のSaturated状態

  対処法:
    - 次回ライブラリ濃度削減: 75 fmol → 50 fmol

原因2: 温度制御異常
  確認方法:
    - MinKNOWで温度警告
    - フローセル温度 >36°C

  対処法:
    - 室温確認、エアコン温度下げ (20-22°C)
    - MinION周辺の通気確保
    - 改善しない: ラン中止、デバイス点検

原因3: 電気的ノイズ
  確認方法:
    - Current histogram: 広範囲のピーク
    - Basecalling品質低下

  対処法:
    - USB接続確認: USB 3.0に直接接続
    - 電源ノイズ除去: UPS使用
    - 周辺電子機器の電源OFF
```

### 9.4 Basecalling遅延問題

**症状: Basecalling queue >5,000 reads**

```yaml
原因1: GPU性能不足
  確認方法:
    - GPU utilization 100%
    - Samples/sec <1,000,000

  対処法:
    a. GPU設定最適化:
       - gpu_runners減少: 8 → 4
       - chunks_per_runner減少: 1024 → 512

    b. オフラインbasecallingに切り替え:
       - MinKNOWでbasecalling無効化
       - FAST5保存を有効化
       - ラン終了後にDorado使用

原因2: ストレージI/O遅延
  確認方法:
    - ディスク使用率 100%
    - 書き込み速度 <100 MB/s

  対処法:
    - 出力先をSSDに変更
    - RAIDストライピング有効化
    - 他プロセスの停止
```

### 9.5 データ品質問題

**症状: Mean Q score <Q15 (24時間後)**

```yaml
原因1: フローセル劣化
  確認方法:
    - フローセル有効期限切れ間近
    - 保管条件不良 (4°C以外で保管)

  対処法:
    - 新しいフローセル使用
    - 保管条件徹底 (4°C、遮光)

原因2: ライブラリ品質不良
  確認方法:
    - TapeStationでアダプターダイマー検出
    - Qubit/TapeStation濃度乖離大

  対処法 (次回改善):
    - AMPure XP精製の徹底
    - アダプターライゲーション条件最適化
    - SFB洗浄回数増加 (1回→2回)

原因3: DNA損傷
  確認方法:
    - Read N50低下 + 品質低下の併発
    - 特定配列領域で品質急降下

  対処法 (次回改善):
    - サンプル保管条件見直し (-80°C徹底)
    - 凍結融解回数削減
    - 抽出時のシェアリング最小化
```

----

## 10. ALCOA+準拠記録

### 10.1 ローディング記録シート

**フローセルローディング記録フォーマット**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MinION フローセルローディング記録
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

実施日時: YYYY-MM-DD HH:MM (ISO 8601形式)
実施者: [フルネーム] _____________ 署名: _______

【1. フローセル情報】
  フローセルID: FAU/FAT番号 __________________
  ロット番号: _________________________________
  有効期限: YYYY-MM-DD
  Platform QC結果:
    アクティブポア数: _______ pores
    判定: □ Excellent  □ Good  □ Acceptable

【2. ライブラリ情報】
  ライブラリID: _______________________________
  サンプルID: __________________________________
  ライブラリタイプ: □ DNA  □ RNA (cDNA)

  濃度測定結果:
    Qubit濃度: _______ ng/μL
    TapeStation N50: _______ bp
    モル濃度: _______ fmol/μL (計算値)

  ローディング量:
    目標fmol: _______ fmol
    使用体積: _______ μL

【3. Loading Mix調製】
  調製日時: YYYY-MM-DD HH:MM

  構成 (チェック):
    □ Sequencing Buffer II: 37.5 μL
    □ ライブラリ: _____ μL
    □ Loading Beads II: 25.5 μL (ボルテックス確認済)
    □ 総量: 75 μL
    □ 気泡なし確認

【4. プライミング実施】
  Flush Buffer調製:
    □ FLT 30 μL + Water 970 μL = 1,000 μL
    調製日時: YYYY-MM-DD HH:MM

  第1回プライミング:
    実施時刻: HH:MM
    □ 吸引 30 μL完了
    □ Flush Buffer 800 μL注入完了
    □ 5分間待機完了

  第2回プライミング (ローディング直前):
    実施時刻: HH:MM
    □ 吸引 200 μL完了
    □ Flush Buffer 200 μL注入完了

【5. ライブラリローディング】
  実施時刻: HH:MM

  手順チェック:
    □ Loading Mix再懸濁 (ピペッティング5回)
    □ プライミングポートから200 μL吸引
    □ Loading Mix 75 μL注入 (30秒以上かけて)
    □ プライミングポートキャップ閉鎖
    □ MinIONカバー閉鎖
    □ LED緑色確認

【6. MinKNOW設定】
  Run name: ____________________________________
  Sequencing kit: SQK-LSK114
  Basecalling mode: □ Duplex (Super accuracy)
  Run duration: _______ hours

  Barcode使用:
    □ なし
    □ あり (Kit: ______________, Barcodes: ____)

【7. シーケンス開始】
  開始時刻: YYYY-MM-DD HH:MM
  予定終了時刻: YYYY-MM-DD HH:MM

  初期メトリクス (開始15分後):
    Active channels: _______
    Reads/sec: _______
    判定: □ 正常  □ 要監視  □ 異常

【8. 環境条件】
  室温: _______ °C
  湿度: _______ %
  MinKNOW version: _______

【9. 特記事項】
  ________________________________________________
  ________________________________________________
  ________________________________________________

【10. 確認・承認】
  記録者: [フルネーム] _____________ 署名: _______
  日時: YYYY-MM-DD HH:MM

  確認者: [フルネーム] _____________ 署名: _______
  日時: YYYY-MM-DD HH:MM

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### 10.2 シーケンスラン監視記録

**ラン中モニタリングシート**

```yaml
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MinION シーケンスラン監視記録
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Run name: ________________________________________
開始日時: YYYY-MM-DD HH:MM
監視者: [フルネーム] _____________

【モニタリング記録】

┌─────────────────────────────────────────────┐
│ 6時間後チェック                                  │
├─────────────────────────────────────────────┤
  日時: YYYY-MM-DD HH:MM
  監視者署名: _______

  メトリクス:
    Active channels: _______
    Total reads: _______
    Total yield: _______ Gb
    Read N50: _______ kb
    Reads/sec: _______
    Mean Q score: Q_______
    Sequencing duty time: _______ %

  判定:
    □ 正常範囲
    □ 要監視 (理由: _________________________)
    □ 異常 (対応: ___________________________)

  特記事項: ______________________________________

└─────────────────────────────────────────────┘

┌─────────────────────────────────────────────┐
│ 24時間後チェック                                 │
├─────────────────────────────────────────────┤
  日時: YYYY-MM-DD HH:MM
  監視者署名: _______

  メトリクス:
    Active channels: _______
    Total reads: _______
    Total yield: _______ Gb
    Read N50: _______ kb
    Duplex rate: _______ %
    Q30+ reads: _______ %
    Mean Q score: Q_______
    Sequencing duty time: _______ %

  判定:
    □ 正常範囲
    □ 要監視 (理由: _________________________)
    □ 異常 (対応: ___________________________)

  早期終了判断:
    □ 継続 (理由: ___________________________)
    □ 終了検討 (理由: ______________________)

  特記事項: ______________________________________

└─────────────────────────────────────────────┘

┌─────────────────────────────────────────────┐
│ 48時間後/ラン終了時チェック                      │
├─────────────────────────────────────────────┤
  日時: YYYY-MM-DD HH:MM
  監視者署名: _______

  最終メトリクス:
    Active channels (終了時): _______
    Total reads: _______
    Total yield: _______ Gb
    Read N50: _______ kb
    Duplex rate: _______ %
    Duplex yield: _______ Gb
    Q30+ reads: _______ %
    Mean Q score: Q_______

  PMDA病原体検出QC判定:
    □ 合格 (Total yield ≥15 Gb, Q30 ≥5 Gb, N50 ≥10 kb)
    □ 条件付き合格 (理由: ___________________)
    □ 不合格 (理由: _________________________)

  データバックアップ:
    □ FASTQ files → NAS
    □ MinKNOW report → ドキュメントサーバー
    □ Run logs → バックアップサーバー

  次のステップ:
    □ バイオインフォマティクス解析へ進む
    □ シーケンス再実施 (理由: _______________)

  特記事項: ______________________________________
  ________________________________________________

└─────────────────────────────────────────────┘

【承認】
  最終確認者: [フルネーム] _____________ 署名: _____
  日時: YYYY-MM-DD HH:MM

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### 10.3 電子記録の管理

**MinKNOWデータの長期保存**

```yaml
保存対象ファイル:

  1. FASTQ files:
     場所: /var/lib/minknow/data/[run_name]/fastq_pass/
     保存期間: 5年以上
     バックアップ: 3箇所 (本番NAS + バックアップNAS + クラウド)

  2. MinKNOW reports:
     - report_[run_name].md
     - report_[run_name].pdf
     - report_[run_name].json
     保存期間: 5年以上

  3. Run logs:
     - duty_time_[run_name].csv
     - throughput_[run_name].csv
     - sequencing_summary_[run_name].txt
     保存期間: 5年以上

  4. Configuration files:
     - run_[run_name]_config.toml
     - user_messages_[run_name].log
     保存期間: 5年以上

バックアップ自動化スクリプト例:

#!/bin/bash
# MinION_data_backup.sh

RUN_NAME=$1
SOURCE_DIR="/var/lib/minknow/data/${RUN_NAME}"
BACKUP_NAS="/mnt/backup_nas/minion_data/${RUN_NAME}"
CLOUD_BACKUP="s3://minion-backup/${RUN_NAME}"

# NASへのバックアップ
echo "Backing up to NAS..."
rsync -avz --progress ${SOURCE_DIR}/ ${BACKUP_NAS}/

# クラウドへのバックアップ
echo "Backing up to Cloud..."
aws s3 sync ${SOURCE_DIR}/ ${CLOUD_BACKUP}/

# チェックサムファイル生成
echo "Generating checksums..."
cd ${SOURCE_DIR}
find . -type f -exec md5sum {} \; > ${RUN_NAME}_checksums.md5

# ログ記録
echo "[$(date)] Backup completed for ${RUN_NAME}" >> \
  /var/log/minion_backup.log

echo "Backup completed successfully."
```

**データ完全性検証**

```yaml
定期検証スケジュール:
  - バックアップ直後: 必須
  - 3ヶ月ごと: 推奨
  - 解析使用前: 必須

検証方法:
  1. ファイル数確認:
     find /path/to/backup -type f | wc -l

  2. チェックサム検証:
     md5sum -c [run_name]_checksums.md5

  3. FASTQファイル完全性:
     zcat *.fastq.gz | head -100000 | \
       awk '{if(NR%4==0) print length($0)}' | \
       sort | uniq -c
     → 品質行の長さが配列行と一致するか確認

  4. MinKNOWレポート再生成可能性:
     MinKNOW CLIでレポート再生成テスト
```

**ALCOA+原則の遵守確認**

```yaml
ALCOA+チェックリスト:

  Attributable (帰属性):
    □ 全記録にフルネーム + 署名
    □ 電子ログにユーザーID記録
    □ タイムスタンプ自動記録

  Legible (判読性):
    □ 手書き記録は明瞭
    □ 電子データは標準フォーマット (FASTQ, CSV)
    □ 略語は定義済み

  Contemporaneous (同時性):
    □ 記録は実施と同時
    □ タイムスタンプ改ざん不可
    □ 遅延記録は理由記載

  Original (原本性):
    □ 生データ保存 (FASTQ)
    □ 修正履歴追跡可能
    □ バックアップと原本の区別明確

  Accurate (正確性):
    □ 測定値は機器から直接転記
    □ 計算値は式と共に記録
    □ エラーは訂正線 + 署名で修正

  Complete (完全性):
    □ 全ステップ記録
    □ 逸脱事項も記録
    □ 欠損データなし

  Consistent (一貫性):
    □ フォーマット統一
    □ 命名規則遵守
    □ 単位表記統一

  Enduring (耐久性):
    □ 5年以上保存可能な媒体
    □ 3箇所以上にバックアップ
    □ 定期的完全性検証

  Available (利用可能性):
    □ 検索可能なインデックス
    □ アクセス権限管理
    □ 監査時に即座提供可能
```

----

## まとめ

本章では、MinIONフローセルのQC、プライミング、ライブラリローディング、シーケンス実行、リアルタイムモニタリング、トラブルシューティング、ALCOA+準拠記録の全プロセスを詳述しました。

**重要ポイント:**

1. **Platform QCの徹底**
   - アクティブポア数 ≥1,200 (最低800)
   - ポア分布の均一性確認
   - 不良フローセルの早期検出

2. **プライミングの正確性**
   - Flush Buffer調製の厳密性
   - 気泡混入の徹底回避
   - 注入速度のコントロール

3. **ライブラリローディング計算**
   - 正確なfmol計算
   - Loading Beads IIの完全懸濁
   - 調製から5分以内のローディング

4. **Duplexベースコーリング設定**
   - Super accuracy (Duplex) mode必須
   - GPU性能の確保
   - Q30以上のデータ取得

5. **リアルタイムQCモニタリング**
   - 6時間、24時間、48時間での定期チェック
   - 異常の早期検出と対処
   - 目標データ量達成の確認

6. **ALCOA+準拠記録**
   - 全ステップの詳細記録
   - 電子データの長期保存
   - データ完全性の定期検証

**次章予告:**

第6章では、MinION用RNAライブラリ調製プロトコルを詳述します。RNA病原体検出のためのcDNA合成、ストランド特異的ライブラリ調製、RNA特異的QC基準について解説します。

----

**文書情報**
- 作成日: 2025-03-08
- バージョン: 1.0
- 作成者: MinIONメタゲノム解析プロトコル開発チーム
- 承認者: [承認者名]
- 次回改訂予定: 2025-09-08
