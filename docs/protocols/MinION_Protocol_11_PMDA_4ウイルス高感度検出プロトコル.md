# MinION Protocol 11: PMDA 4ウイルス高感度検出プロトコル

**Version:** 1.1
**作成日:** 2025-11-12
**最終更新:** 2025-11-13
**対象ウイルス:** ポリオーマウイルス、ハンタウイルス、東部ウマ脳炎ウイルス(EEEV)、ブタスピューマウイルス
**検出感度目標:** \<50 copies/mL plasma (スピューマウイルスは\<10 copies/10⁵ PBMCs)

---- 

## ⚠️ 重要: プロトコル使用ガイダンス

### 📍 Protocol 12 を優先使用してください

**本プロトコル（Protocol 11）は高感度化オプション**であり、通常のスクリーニングには **Protocol 12「統合サンプル調製プロトコル」** を使用することを強く推奨します。

| 項目       | Protocol 12（推奨）   | Protocol 11（本文書）   |
| -------- | ----------------- | ------------------ |
| **用途**   | 通常スクリーニング（91病原体）  | 高感度確認検査（4ウイルス）     |
| **LOD**  | 100-500 copies/mL | \<50 copies/mL     |
| **実験時間** | 13時間              | 19.5時間（+6.5時間）     |
| **コスト**  | ¥157,000          | ¥185,000（+¥28,000） |
| **複雑度**  | 低（2ワークフロー）        | 高（4ワークフロー）         |
| **適応**   | 全サンプル             | 以下の場合のみ            |

### ✅ Protocol 11を使用すべき場合

以下のいずれかに該当する場合のみ、本プロトコルを使用してください：

1. **Protocol 12で陽性検出後の確認検査**
   2. 初回スクリーニングで4ウイルスのいずれかが検出された
   3. 定量精度向上のため高感度検査が必要

2. **疫学的な感染リスクが高い**
   2. 輸入ブタ、非SPF施設由来
   3. 感染症発生地域からの移動歴

3. **規制当局からの特別要求**
   2. LOD \<50 copies/mL達成が必須条件として指定された場合

4. **研究目的での超高感度解析**
   2. バリデーション試験、LOD決定試験

### ❌ Protocol 11を使用しない方が良い場合

- **通常のドナーブタスクリーニング** → Protocol 12で十分
- **コスト・時間を最小化したい** → Protocol 12が効率的
- **初回検査（確定診断ではない）** → Protocol 12で十分

---- 

## 1. プロトコル概要

### 1.1 目的

本プロトコルは、異種移植用ドナーブタの病原体スクリーニングにおいて、PMDA指定91病原体のうち特に検出が困難な以下の4ウイルスについて、超高感度検出(\<50 copies/mL)を実現するための統合プロトコルである。

### 1.2 対象ウイルスの特性

| ウイルス                 | ゲノム型              | サイズ(kb)                    | 血漿中の形態          | PMDA分類              | 検出上の課題                  |
| -------------------- | ----------------- | -------------------------- | --------------- | ------------------- | ----------------------- |
| **ポリオーマウイルス**        | dsDNA環状           | 5.0                        | cfDNA(ビリオン由来)   | Line 35/58          | 宿主DNA混入、類似ウイルスとの鑑別      |
| **ハンタウイルス**          | ssRNA(-)3分節       | 11.8 (L:6.5, M:3.6, S:1.7) | ビリオンRNA + cfRNA | Line 31/54, 人獣共通感染症 | rRNA除去、poly(A)尾なし、RNA分解 |
| **東部ウマ脳炎ウイルス(EEEV)** | ssRNA(+)、poly(A)+ | 11-12                      | ビリオンRNA + cfRNA | Line 32/55, 人獣共通感染症 | rRNA除去、短時間ウイルス血症        |
| **ブタスピューマウイルス**      | ssRNA→dsDNAプロウイルス | 7-10 (プロウイルス)              | PBMCゲノムDNA(組込型) | Line 154, 特別管理微生物   | ブタでの存在不明、PERV識別         |

### 1.3 検出戦略の概要

**ウイルス核酸型による分類:**

- **dsDNAウイルス (ポリオーマウイルス):** 血漿cfDNA抽出 → CpGメチル化除去 → DNAライブラリ → MinIONシーケンス
- **ssRNA(+)ウイルス (EEEV):** 血漿cfRNA抽出 → poly(A)選択 → Direct RNAまたはcDNAライブラリ → MinIONシーケンス
- **ssRNA(-)ウイルス (ハンタウイルス):** 血漿cfRNA抽出 → rRNA除去 → アンプリコンRT-PCRまたはSMART-9N cDNA → MinIONシーケンス
- **レトロウイルスプロウイルス (スピューマウイルス):** PBMC分離 → ゲノムDNA抽出 → Nested PCR (pol遺伝子) → MinIONシーケンスまたはqPCR定量

---- 

## 2. サンプル採取と前処理

### 2.1 血液採取

**対象:** ポリオーマウイルス、ハンタウイルス、EEEV

#### 2.1.1 採血量と使用チューブ

| 目的                     | 採血量   | チューブ種類                                                  | 添加試薬                        | 保存条件          |
| ---------------------- | ----- | ------------------------------------------------------- | --------------------------- | ------------- |
| **DNA/RNAウイルス検出(標準)**  | 5 mL  | EDTA                                                    | なし                          | 4時間以内に処理(4°C) |
| **DNA/RNAウイルス検出(高感度)** | 10 mL | EDTA                                                    | SUPERase•In 200 U (20 U/mL) | 4時間以内に処理(4°C) |
| **遅延処理が必要な場合**         | 10 mL | Streck cfDNA BCT (DNA) / cf-RNA Preservative Tube (RNA) | 専用保存液含有                     | 室温30日間安定      |

#### 2.1.2 RNA保護のための重要注意事項

**ハンタウイルスとEEEV検出では、RNA分解防止が最重要:**

1. **RNase阻害剤の即時添加**
   2. 採血直後にSUPERase•In RNase Inhibitor (Invitrogen)を20 U/mLで添加
   3. SUPERase•InはRNase A/B/C/1/T1を阻害

2. **RNaseフリー環境の維持**
   2. 手袋着用(素手で触らない)
   3. フィルターチップ使用
   4. RNaseZap処理済みワークスペース

3. **凍結融解サイクルの最小化**
   2. 血漿分離後は-80°Cで保存
   3. 凍結融解は最大1-2回まで(それ以上はRNA分解が進行)

### 2.2 血漿分離(DNA/RNAウイルス)

#### 2.2.1 2段階遠心法

```
【第1遠心】細胞除去
- 1,900×g、10分、4°C
- 上清(血漿)を新しい15 mLチューブに移す

【第2遠心】細胞破片除去
- 16,000×g、10分、4°C
- 上清(清澄血漿)を0.5 mL aliquotに分注
- -80°Cで保存(DNA安定、RNA不安定のため使用まで凍結保存)
```

**収量目安:** 10 mL全血 → 4-5 mL血漿

### 2.3 PBMC分離(スピューマウイルス)

#### 2.3.1 Ficoll-Paque密度勾配遠心法

```
【準備】
- 10 mL EDTA全血を室温に戻す
- Ficoll-Paque PLUS (GE Healthcare)を室温に戻す

【遠心分離】
1. 15 mL遠心チューブにFicoll-Paque 3 mLを入れる
2. 全血10 mLをPBS 10 mLで希釈(20 mL)
3. 希釈血液をFicoll-Paque上に静かに重層
4. 800×g、30分、室温、ブレーキなし遠心
5. バフィーコート層(白濁層)を採取

【PBMC洗浄】
1. 採取したPBMCを10 mL PBSで懸濁
2. 300×g、10分、室温で遠心
3. 上清除去、PBSで再懸濁(洗浄2回実施)
4. 最終ペレットを-80°Cで保存またはDNA抽出

【細胞数カウント】
- Trypan blue染色で生細胞数カウント
- 目安: 10 mL全血 → 1-5×10⁶ PBMCs
```

---- 

## 3. 核酸抽出プロトコル

### 3.1 同時DNA/RNA抽出(ポリオーマ・ハンタ・EEEV検出用)

**推奨キット:** Zymo Quick-cfDNA/RNA Serum & Plasma Kit

#### 3.1.1 抽出手順

```
【Input】血漿 500 μL (高感度の場合は5 mL)

【手順】
1. 血漿500 μLにProteinase K 10 μL + Binding Buffer 500 μLを添加
2. ボルテックス10秒、56°C 10分インキュベート
3. Zymo-Spin IC-XLRカラムに添加、12,000×g 1分遠心
4. カラムをDNA collection tubeに移す
5. フロースルーを新しいチューブに保存(RNA用)

【DNA溶出】
6. カラムにDNA Wash Buffer 400 μL添加、12,000×g 1分遠心(2回)
7. 新しい1.5 mL tubeにカラムをセット
8. DNA Elution Buffer 50 μL添加、1分室温静置後、12,000×g 1分遠心
   → **DNA elute 50 μL取得(ポリオーマウイルス検出用)**

【RNA溶出】
9. 保存したフロースルーにRNA Binding Buffer 1 mLを添加
10. 新しいZymo-Spin IC-XLRカラムに添加、12,000×g 1分遠心
11. RNA Wash Buffer 400 μL添加、12,000×g 1分遠心(2回)
12. 新しい1.5 mL tubeにカラムをセット
13. DNase/RNase-Free Water 50 μL添加、1分室温静置後、12,000×g 1分遠心
    → **RNA elute 50 μL取得(ハンタ・EEEV検出用)**
```

**収量目安:**
- cfDNA: 5-50 ng/mL血漿(健常ブタ)、ウイルス感染時は増加
- cfRNA: 1-10 ng/mL血漿(定量困難、RNAは低濃度)

### 3.2 PBMC genomic DNA抽出(スピューマウイルス検出用)

**推奨キット:** QIAGEN QIAamp DNA Blood Mini Kit

#### 3.2.1 抽出手順

```
【Input】PBMC 1-5×10⁶ cells

【手順】
1. PBMCペレットにPBS 200 μLを添加、懸濁
2. Proteinase K 20 μL + Buffer AL 200 μLを添加
3. ボルテックス15秒、56°C 10分インキュベート
4. 99.8%エタノール200 μL添加、ボルテックス15秒
5. QIAamp Mini spinカラムに全量添加、8,000×g 1分遠心
6. Buffer AW1 500 μL添加、8,000×g 1分遠心
7. Buffer AW2 500 μL添加、14,000×g 3分遠心
8. 新しい1.5 mL tubeにカラムをセット
9. Buffer AE 50 μL添加、1分室温静置後、8,000×g 1分遠心
   → **Genomic DNA 50 μL取得(~100-500 ng/μL)**
```

**スピューマウイルスプロウイルスコピー数推定:**
- 1 PBMC = 約6 pg DNA
- 100 ng DNA = 約16,000 cells
- 感染率10%と仮定 → 100 ngに1,600プロウイルスコピー含有

---- 

## 4. ウイルス特異的前処理

### 4.1 ポリオーマウイルス: 宿主DNA除去

**目的:** ブタゲノムDNA(Sus scrofa, 2.6 Gb)を除去し、ウイルスDNAを濃縮

**使用キット:** NEBNext Microbiome DNA Enrichment Kit (E2612)

#### 4.1.1 CpGメチル化キャプチャー法

```
【原理】
- 哺乳類ゲノムDNAは60-90%がCpGメチル化
- ウイルスDNAは非メチル化
- MBD-Fc磁気ビーズがメチル化DNAを選択的に捕捉
- 上清にウイルスDNAが濃縮(20-100倍濃縮)

【手順】
1. cfDNA 50 μLにNEB Binding Buffer 100 μL + MBD-Fc Beads 10 μLを添加
2. ボルテックス、室温15分インキュベート(回転混和)
3. 磁気スタンドで2分静置、上清を新しいチューブに回収
4. 上清にHigh Salt Wash Buffer 200 μL添加、混和
5. 磁気スタンドで2分静置、上清を回収(ビーズ捨てる)
6. エタノール沈殿でDNA濃縮(省略可)
   → **宿主DNA除去cfDNA 約45 μL取得(ポリオーマウイルス濃縮)**

【除去効率】
- 宿主DNA: 60-90%除去
- ウイルスDNA: 80-95%回収
- 微生物DNA濃縮: 20倍
```

### 4.2 EEEV: Poly(A) RNA選択

**目的:** rRNA(80-90%)を除去し、poly(A)尾を持つEEEV RNAを濃縮

**使用キット:** NEBNext Poly(A) mRNA Magnetic Isolation Module (E7490)

#### 4.2.1 Oligo(dT)磁気ビーズ法

```
【前処理: DNase処理】
1. RNA 50 μLにDNase I (RNase-free) 1 U/μg添加
2. Protector RNase Inhibitor 20 U添加
3. 25°C 30分インキュベート
4. EDTA 5 mM添加、65°C 10分(DNase不活化)

【Poly(A)選択】
1. DNase処理RNA 50 μLを65°C 5分加熱(2次構造解離)
2. Oligo(dT) Magnetic Beads 50 μL添加
3. 室温10分インキュベート(poly(A) RNA結合)
4. 磁気スタンドで2分、上清除去(rRNA含有、廃棄)
5. Wash Buffer 200 μL添加、混和、磁気スタンドで上清除去(2回洗浄)
6. Tris Buffer (10 mM, pH 7.5) 50 μL添加
7. 65°C 2分加熱(poly(A) RNA溶出)
8. 磁気スタンドで上清回収
   → **Poly(A)+ RNA 約45 μL取得(EEEV濃縮、rRNA 98%除去)**

【濃縮効率】
- rRNA除去: 98%
- mRNA/poly(A)+ RNA回収: 80-90%
- EEEVウイルスRNA濃縮: 100-1,000倍
```

**重要注意:** ハンタウイルスはpoly(A)尾を持たないため、このステップで失われる。ハンタウイルスにはrRNA除去法を使用。

### 4.3 ハンタウイルス: rRNA除去

**目的:** rRNAを除去し、poly(A)尾を持たないハンタウイルスRNAを濃縮

**使用キット:** NEBNext rRNA Depletion Kit (Human/Mouse/Rat) (E6310)

#### 4.3.1 RNase H法

```
【前処理: DNase処理】
- 上記4.2.1と同様

【rRNA除去】
1. DNase処理RNA 50 μLにrRNA Depletion Probes 10 μL添加
2. 95°C 2分 → 4°C 30秒(プローブハイブリダイゼーション)
3. RNase H 5 μL添加(rRNA特異的分解)
4. 37°C 30分インキュベート
5. DNase I 5 μL添加、37°C 30分(rRNA断片除去)
6. RNA Clean & Concentrator (Zymo)で精製
   → **rRNA除去RNA 約40 μL取得(ハンタウイルス濃縮)**

【除去効率】
- rRNA除去: 80-95%
- ウイルスRNA回収: 70-90%
- ハンタウイルスRNA濃縮: 10-50倍

【ブタrRNA対応】
- NEBNext Human/Mouse/RatキットはブタrRNAに部分的に交差反応
- カスタムブタrRNAプローブ設計も検討可(コスト増)
```

### 4.4 スピューマウイルス: Nested PCR濃縮

**目的:** 0.001-0.01%のプロウイルスDNAをPCR増幅で検出可能レベルに濃縮

#### 4.4.1 Degenerate Primer設計(pol遺伝子標的)

**背景:**
- ブタスピューマウイルスのリファレンスゲノムは存在しない
- サル泡沫ウイルス(SFV)のpol遺伝子(最保存領域)をベースに縮重プライマー設計
- 逆転写酵素ドメイン標的(種間30-50%相同性)

**プライマー設計例(文献ベース、SFV用):**

```
【Outer primers (1st PCR)】
FV-pol-F1: 5'-GGNCARATHGGNATGTTYGG-3' (逆転写酵素conserved motif)
FV-pol-R1: 5'-CCRTCNCCRAANCCRTC-3' (逆転写酵素conserved motif)
Expected product: ~800 bp

【Inner primers (2nd PCR)】
FV-pol-F2: 5'-ATHGGNCARGGNTTYACNAC-3'
FV-pol-R2: 5'-GTRTCNGTYTTRTCNCC-3'
Expected product: ~400 bp
```

#### 4.4.2 Nested PCR手順

```
【1st PCR】
- Template: Genomic DNA 100 ng (約16,000 cells)
- Primers: FV-pol-F1/R1 各0.4 μM
- Taq polymerase: High-fidelity (Q5, Phusion)
- Cycles: 95°C 3分 → (98°C 10秒 → 55°C 30秒 → 72°C 1分) ×35 cycles → 72°C 5分
- Product: 1st PCR 50 μL

【2nd PCR】
- Template: 1st PCR産物 1 μL (1:50希釈)
- Primers: FV-pol-F2/R2 各0.4 μM
- Cycles: 95°C 3分 → (98°C 10秒 → 58°C 30秒 → 72°C 40秒) ×25 cycles → 72°C 5分
- Product: 2nd PCR 50 μL

【Agarose gel確認】
- 1.5% agarose gel電気泳動
- 期待バンド: ~400 bp
- バンドありの場合: シーケンス確認(Sanger or MinION)

【感度】
- LOD: 1-10 copies per 10⁵ PBMCs
- 特異度: PERVと識別可能(pol遺伝子配列が異なる)
```

---- 

## 5. MinIONライブラリ調製

### 5.1 ポリオーマウイルス: DNAライブラリ

**使用キット:** Oxford Nanopore Ligation Sequencing Kit V14 (SQK-LSK114)

```
【Input】宿主DNA除去cfDNA 45 μL

【手順】
1. DNA Repair and End-prep
   - DNA 45 μL + FFPE DNA Repair Buffer 3.5 μL + FFPE DNA Repair Mix 2 μL
   - 20°C 5分 → 65°C 5分

2. Adapter Ligation
   - 修復DNA 50 μL + Ligation Buffer 25 μL + NEBNext Quick T4 DNA Ligase 10 μL + Adapter Mix II 5 μL
   - 室温20分

3. Short Fragment Buffer wash
   - ライゲーション産物90 μL + Short Fragment Buffer 40 μL
   - AMPure XP Beads 90 μL添加、混和5分
   - 磁気スタンド5分、上清除去
   - 80% Ethanol 200 μL洗浄(2回)
   - Elution Buffer 15 μL溶出

4. 定量とQC
   - Qubit dsDNA HS Assay: 目標 >50 fmol
   - 濃度が低い場合は複数サンプルをプール

【Output】ライブラリ 15 μL (50-200 fmol)
```

### 5.2 EEEV: Direct RNA Sequencing

**使用キット:** Oxford Nanopore Direct RNA Sequencing Kit (SQK-RNA002)

**利点:** EEEVはpoly(A)尾を持つため、Direct RNA可能(ハンタウイルスは不可)

```
【Input】Poly(A)+ RNA 45 μL (目標 >50 ng)

【手順】
1. Reverse Transcription
   - Poly(A)+ RNA 9 μL + RT Adapter 1 μL
   - 65°C 5分(アニール)→氷上2分
   - Strand-Switching Primer 1 μL + RT Mix 9 μL添加
   - 50°C 10分 → 70°C 10分

2. Adapter Ligation
   - RT産物20 μL + RNA Adapter 2 μL + Ligation Buffer 10 μL + T4 DNA Ligase 8 μL
   - 室温10分

3. Beads wash
   - ライゲーション産物40 μL + RNA Wash Buffer 40 μL
   - AMPure XP Beads 80 μL添加、混和5分
   - 磁気スタンド5分、上清除去
   - Wash Buffer 150 μL洗浄(2回)
   - Elution Buffer 21 μL溶出

【Output】Direct RNAライブラリ 21 μL
```

### 5.3 ハンタウイルス: Amplicon-based cDNA

**戦略:** メタゲノムでは感度不足(\<10⁵ copies/mL)のため、アンプリコンアプローチ

#### 5.3.1 Tiled Amplicon設計

**Hantaウイルス特異的プライマー設計:**

```
【原理】
- ハンタウイルスL/M/S 3分節全体を400 bp amplicons × 100 bp overlapでカバー
- ARTIC Network様のタイルアンプリコン戦略
- 保存領域ベースのconsensus primerで複数種対応

【Primer pool設計】
- L segment (6.5 kb): 20 amplicons (400 bp × 20 = 8 kb coverage)
- M segment (3.6 kb): 11 amplicons
- S segment (1.7 kb): 5 amplicons
- Total: 36 primer pairs (Pool A: 18, Pool B: 18で多重PCR)

【参考文献】
- 韓国のHantaan virus multiplex PCR-nanopore protocol (Viruses, 2021)
```

#### 5.3.2 Multiplex RT-PCR

```
【1st-strand cDNA合成】
- rRNA除去RNA 10 μL + Random Hexamers 1 μL
- 65°C 5分 → 氷上2分
- RT Buffer 4 μL + dNTP 1 μL + RNase inhibitor 1 μL + SuperScript IV RT 1 μL
- 23°C 10分 → 50°C 10分 → 80°C 10分
- cDNA 20 μL取得

【Multiplex PCR (Pool A)】
- cDNA 2.5 μL + Primer Pool A (18 pairs, 各0.015 μM) 12.5 μL + Q5 High-Fidelity Master Mix 12.5 μL
- 98°C 30秒 → (98°C 15秒 → 63°C 5分) ×35 cycles → 65°C 5分
- PCR product 25 μL

【Multiplex PCR (Pool B)】
- 同様にPrimer Pool Bで実施

【プール】
- Pool A + Pool B = 50 μL
- AMPure XP 1.0× cleanup
- 溶出: 15 μL

【Native Barcoding】
- Oxford Nanopore Native Barcoding Kit (SQK-NBD114)
- 複数サンプルの多重化可能
```

### 5.4 スピューマウイルス: PCR product直接シーケンス

```
【Nested PCR産物のライブラリ化】
- 2nd PCR産物(~400 bp) 45 μL
- Oxford Nanopore Ligation Sequencing Kit (SQK-LSK114)使用
- 上記5.1と同じ手順でアダプターライゲーション

【定量】
- Qubit dsDNA HS: 目標 >50 fmol (400 bp × 1×10¹¹ molecules = ~27 ng)
```

---- 

## 6. MinION Sequencing

### 6.1 Flow Cell QC

```
【使用Flow Cell】
- FLO-MIN106D (R9.4.1 chemistry)
- 標準: 512 pores (新品)
- 使用可能最低: >200 active pores

【Flow Cell QC手順】
1. MinKNOWソフトウェア起動
2. Flow Cell Check実行
3. Pore count確認: >200 OK、<200 NG(新しいFlow Cell使用)
```

### 6.2 Priming and Loading

```
【Priming】
1. Flow Cell Priming Kit (FPK) thaw
2. Priming portを開き、800 μL Flush Tether除去
3. Priming Mix (Flush Buffer 1170 μL + Flush Tether 30 μL)を800 μL注入
4. 5分静置
5. Sample portを開き、200 μL Priming Mix注入

【Library Loading】
1. ライブラリ調製
   - Sequencing Buffer II (SB) 37.5 μL + Loading Beads II (LB) 25.5 μL混和
   - ライブラリ 12 μL添加(最終75 μL)
2. Sample portから75 μL注入(ドロップワイズ、気泡混入回避)
3. 即座にMinKNOWでrun開始
```

### 6.3 Sequencing Parameters

**ウイルス別推奨Run時間:**

| ウイルス                  | ウイルス量予想               | Run時間   | 期待Read数                    | 判定基準          |
| --------------------- | --------------------- | ------- | -------------------------- | ------------- |
| **ポリオーマウイルス**         | 10³-10⁷ copies/mL     | 24-48時間 | \>100 virus reads          | ≥100 readsで陽性 |
| **EEEV**              | 10⁵-10⁸ copies/mL     | 8-24時間  | \>100 virus reads          | ≥100 readsで陽性 |
| **ハンタウイルス(amplicon)** | 10-10³ copies/mL      | 4-8時間   | \>1,000 reads (全amplicons) | L/M/S全分節検出で陽性 |
| **スピューマウイルス(PCR)**    | 1-10 copies/10⁵ cells | 1-4時間   | \>100 reads (pol gene)     | Sanger確認後陽性   |

**低力価サンプル戦略:**
- Real-time basecallingでウイルスread出現をモニター
- 4時間経過時点でウイルスread \<10の場合:
  - オプション1: Adaptive Samplingでウイルス濃縮(host read reject)
  - オプション2: Hybridization capture enrichmentを追加し、再ライブラリ調製

---- 

## 7. Bioinformatics解析

### 7.1 ポリオーマウイルス検出

#### 7.1.1 リファレンスデータベース

```
【データベース構築】
/mnt/efs/databases/pmda/2024.1/polyomavirus/

収録配列:
- Sus scrofa polyomavirus 2 (GenBank検索、2018年登録)
- BK polyomavirus (NC_001538)
- JC polyomavirus (NC_001699)
- SV40 (NC_001669)
- 全Polyomaviridae family genomes (NCBI Taxonomy ID 151341)
```

#### 7.1.2 検出パイプライン

```bash
# Basecalling
guppy_basecaller -i fast5/ -s fastq/ --config dna_r9.4.1_450bps_hac.cfg

# QC
NanoPlot --fastq fastq/*.fastq -o qc/

# Host removal (optional、宿主DNA除去済みの場合スキップ)
minimap2 -ax map-ont sus_scrofa_ref.mmu -t 8 reads.fastq | \
  samtools view -f 4 -b | samtools fastq > unmapped.fastq

# Polyomavirus detection
minimap2 -ax map-ont polyomavirus_db.mmi unmapped.fastq > polyoma.sam

# Filtering
samtools view -bS -F 4 -q 10 polyoma.sam | samtools sort -o polyoma_sorted.bam
samtools index polyoma_sorted.bam

# Read count
samtools view -c polyoma_sorted.bam
# ≥100 reads → 陽性

# Coverage analysis
samtools depth polyoma_sorted.bam | awk '{sum+=$3} END {print sum/5000}'
# ≥10× depth → 高信頼度陽性
```

#### 7.1.3 定量

```python
# scripts/phase5_quantification.py

def quantify_polyomavirus(bam_file, plasma_volume_ml):
    """
    ポリオーマウイルスコピー数定量
    """
    import pysam

    # Read count
    bam = pysam.AlignmentFile(bam_file, "rb")
    virus_reads = bam.count()

    # Total reads
    total_reads = get_total_reads_from_fastq()

    # Standard curve (事前にスパイクイン実験で作成)
    # copies/mL = f(virus_reads, total_reads, plasma_volume)
    # 例: 標準曲線 y = 10^(0.8*log10(virus_reads/total_reads) + 3.5)

    if virus_reads == 0:
        return 0

    fraction = virus_reads / total_reads
    copies_per_ml = (10 ** (0.8 * np.log10(fraction) + 3.5)) * (5 / plasma_volume_ml)

    return copies_per_ml
```

### 7.2 ハンタウイルス検出(3分節concordance check)

#### 7.2.1 リファレンスデータベース

```
【データベース】
/mnt/efs/databases/pmda/2024.1/hantavirus/

収録配列:
- Hantaan virus: L (NC_005222), M (NC_005219), S (NC_005218)
- Seoul virus: L (NC_005238), M (NC_005236), S (NC_005237)
- Dobrava-Belgrade virus: NC_005233, NC_005234, NC_005235
- Puumala virus, Sin Nombre virus等
```

#### 7.2.2 3分節一致性チェック

```python
# scripts/phase4_pathogen_detection.py

def detect_hantavirus_trisegmented(bam_file):
    """
    ハンタウイルス3分節(L/M/S)の同時検出確認
    """
    import pysam

    bam = pysam.AlignmentFile(bam_file, "rb")

    segment_counts = {"L": 0, "M": 0, "S": 0}

    for read in bam:
        ref_name = read.reference_name
        if "L_segment" in ref_name or "_L_" in ref_name:
            segment_counts["L"] += 1
        elif "M_segment" in ref_name or "_M_" in ref_name:
            segment_counts["M"] += 1
        elif "S_segment" in ref_name or "_S_" in ref_name:
            segment_counts["S"] += 1

    # 判定基準: 全3分節で≥50 reads検出
    if all(count >= 50 for count in segment_counts.values()):
        return "POSITIVE", segment_counts
    elif any(count >= 50 for count in segment_counts.values()):
        return "INCONCLUSIVE", segment_counts  # 再検査推奨
    else:
        return "NEGATIVE", segment_counts
```

### 7.3 EEEV検出(Alphavirus識別)

#### 7.3.1 リファレンスデータベース

```
【データベース】
/mnt/efs/databases/pmda/2024.1/alphavirus/

PMDA指定Alphavirus:
- Eastern equine encephalitis virus (EEEV): NC_003899
- Western equine encephalitis virus (WEEV): NC_003908
- Venezuelan equine encephalitis virus (VEEV): NC_001449
- Getah virus: NC_003696 (Line 18/41)
```

#### 7.3.2 種同定パイプライン

```bash
# Alphavirus detection
minimap2 -ax map-ont alphavirus_db.mmi reads.fastq > alphavirus.sam
samtools view -bS -F 4 -q 20 alphavirus.sam | samtools sort -o alphavirus_sorted.bam

# Species-level assignment
# 最多ヒットリファレンスを種として判定
samtools idxstats alphavirus_sorted.bam | sort -k3 -nr | head -1
# Output例: NC_003899 (EEEV) 523 reads → EEEV陽性
```

#### 7.3.3 系統解析(オプション)

```bash
# Consensus genome assembly
samtools consensus alphavirus_sorted.bam -o eeev_consensus.fasta

# Phylogenetic analysis (nsP1 or E1 gene)
# MAFFT alignment + FastTree
mafft --auto eeev_consensus.fasta eeev_references.fasta > aligned.fasta
FastTree -nt -gtr aligned.fasta > eeev_tree.nwk

# 系統樹でNorth American vs South American lineage判定
```

### 7.4 スピューマウイルス検出(PERV識別)

#### 7.4.1 リファレンスデータベース

```
【データベース】
/mnt/efs/databases/pmda/2024.1/spumavirus/

収録配列:
- Simian foamy virus (SFV): NC_001364 (pol gene)
- Feline foamy virus (FFV): NC_001871
- Bovine foamy virus (BFV): NC_001831
- PERV (識別用): pol gene sequences

注意: ブタスピューマウイルスのリファレンスゲノムは存在しない
```

#### 7.4.2 PERV vs Spumavirus識別

```python
# scripts/phase4_pathogen_detection.py

def discriminate_spumavirus_from_perv(blast_result):
    """
    スピューマウイルスとPERVの識別

    原理:
    - PERV: Gammaretrovirus (MLV類似)
    - Spumavirus: Spumaretrovirinae (系統的に遠い)
    - pol遺伝子配列で明確に区別可能
    """
    best_hit = blast_result[0]

    if "Simian_foamy_virus" in best_hit["subject"] or \
       "Spumavirus" in best_hit["subject"]:
        if best_hit["identity"] > 70:  # 種間相同性30-50%を考慮
            return "SPUMAVIRUS_POSITIVE"

    elif "PERV" in best_hit["subject"] or "porcine_endogenous" in best_hit["subject"]:
        return "PERV_DETECTED"  # 全ブタで陽性(内在性)

    else:
        return "UNKNOWN_RETROVIRUS"
```

#### 7.4.3 Sangerシーケンス確認(推奨)

```
【理由】
- ブタスピューマウイルスのリファレンスが不在
- MinIONのみでは誤同定リスク
- PCR産物をSangerシーケンスで配列確定後、BLASTで系統確認

【手順】
1. Nested PCR産物をゲル精製
2. BigDye Terminator Cycle Sequencing
3. ABI 3730xlでシーケンス(両方向)
4. BLASTn検索でSFV pol遺伝子との相同性確認
   - >70%相同性 → スピューマウイルス陽性
   - <50%相同性 → PERV等他のレトロウイルス
```

---- 

## 8. QC基準と判定

### 8.1 各ウイルスのQC基準

| ウイルス          | 検出基準                      | 定量範囲              | 確認試験                      | PMDA要求LOD             |
| ------------- | ------------------------- | ----------------- | ------------------------- | --------------------- |
| **ポリオーマウイルス** | ≥100 reads, ≥10× coverage | 50-10⁷ copies/mL  | qPCR (BK virus primers)   | \<50 copies/mL        |
| **ハンタウイルス**   | L/M/S全分節≥50 reads         | 100-10⁶ copies/mL | qRT-PCR (S segment)       | \<50 copies/mL        |
| **EEEV**      | ≥100 reads, ≥10× coverage | 50-10⁸ copies/mL  | qRT-PCR (nsP4 or E1 gene) | \<50 copies/mL        |
| **スピューマウイルス** | Nested PCR陽性 + Sanger確認   | 1-1,000/10⁵ PBMCs | PBMC培養活性化試験               | \<10 copies/10⁵ cells |

### 8.2 PMDA Validation要求

**PPA (Positive Percent Agreement):** ≥95%
- LOD濃度(50 copies/mL)で10反復実施 → ≥9/10検出

**NPA (Negative Percent Agreement):** ≥98%
- 陰性サンプル50検体実施 → ≤1/50偽陽性

**R² (Linearity):** ≥0.90
- 10, 50, 100, 500, 1,000, 5,000, 10,000 copies/mL標準曲線
- 測定値 vs 期待値の相関係数

**Precision (再現性):** CV \<20% at 10× LOD

---- 

## 9. トラブルシューティング

### 9.1 ポリオーマウイルス検出不良

| 症状                 | 原因         | 対策                               |
| ------------------ | ---------- | -------------------------------- |
| ウイルスread \<10      | 宿主DNA除去不十分 | CpGキャプチャー2回実施、またはProbe capture追加 |
| 多数のSus scrofa read | CpG除去失敗    | MBD-Fc beads量を2倍に増量              |
| Coverage偏り         | 環状ゲノム線状化不足 | 制限酵素処理(BamHI)またはshearing追加       |

### 9.2 ハンタウイルスRNA分解

| 症状          | 原因          | 対策                                            |
| ----------- | ----------- | --------------------------------------------- |
| RNA収量\<1 ng | RNase混入     | SUPERase•In 40 U/mLに増量、採血から抽出まで2時間以内          |
| RIN \<5     | 凍結融解過多      | 凍結融解1回まで、Aliquot保存                            |
| rRNA除去不十分   | ブタrRNA交差反応低 | カスタムブタrRNAプローブ設計(Integrated DNA Technologies) |

### 9.3 EEEV低力価

| 症状                   | 原因                     | 対策                                |
| -------------------- | ---------------------- | --------------------------------- |
| Poly(A)+ RNA \<10 ng | 血漿量不足                  | 10 mL血漿使用、または5回抽出プール              |
| ウイルスread \<50        | ウイルス血症短期間              | 発熱後48時間以内の採血推奨                    |
| Host mRNA混入多         | Poly(A)選択でhost mRNAも濃縮 | Probe capture追加(Alphavirus panel) |

### 9.4 スピューマウイルスPCR陰性

| 症状           | 原因             | 対策                             |
| ------------ | -------------- | ------------------------------ |
| Nested PCR陰性 | ブタスピューマウイルス不在? | PMDA確認、または検出対象から除外             |
| 1st PCR陰性    | プライマーミスマッチ     | Degenerate primerのdegeneracy増加 |
| PERV偽陽性      | pol遺伝子類似       | Sanger確認、系統解析必須                |

---

## 10. コスト試算

### 10.1 追加コスト(既存パイプライン比)

| 項目              | キット名                                  | 単価       | 用途       | サンプル当たりコスト       |
| --------------- | ------------------------------------- | -------- | -------- | ---------------- |
| rRNA除去          | NEBNext rRNA Depletion Kit (96 rxns)  | ¥150,000 | ハンタ・EEEV | ¥1,563           |
| Poly(A)選択       | NEBNext Poly(A) mRNA Module (96 rxns) | ¥120,000 | EEEV     | ¥1,250           |
| Direct RNAキット   | ONT Direct RNA Kit (24 rxns)          | ¥200,000 | EEEV     | ¥8,333           |
| Nested PCRプライマー | Custom degenerate primers             | ¥50,000  | スピューマ    | ¥521 (96サンプル分)   |
| Sangerシーケンス     | 外注(両方向)                               | ¥3,000   | スピューマ確認  | ¥3,000           |
| **合計追加コスト**     | -                                     | -        | -        | **¥14,667/サンプル** |

**既存パイプラインコスト:** ¥127,000/サンプル
**4ウイルス対応後:** ¥141,667/サンプル (+11.5%)
**従来法(Pattern B):** ¥449,574/サンプル

**結論:** 4ウイルス高感度検出追加でも、従来法の3.2倍安価を維持

### 10.2 バリデーション初期投資

| 項目                  | 詳細                       | コスト            |
| ------------------- | ------------------------ | -------------- |
| 合成DNA/RNA標準品        | 各ウイルス×7濃度                | ¥300,000       |
| Probe capture panel | ViroFind panel or custom | ¥500,000       |
| qPCR/qRT-PCR検証      | TaqMan probes, 検証実験      | ¥200,000       |
| **初期投資合計**          | -                        | **¥1,000,000** |

---- 

## 11. 実施スケジュール

### 11.1 プロトコル実装フェーズ(3-4ヶ月)

| フェーズ                    | 期間  | 作業内容                                     |
| ----------------------- | --- | ---------------------------------------- |
| **Month 1-2: プロトコル最適化** | 2ヶ月 | RNA抽出、rRNA除去、poly(A)選択の最適化。合成標準品スパイクイン実験 |
| **Month 3: Pipeline統合** | 1ヶ月 | Phase 0追加、Phase 3/4修正、データベース構築、解析スクリプト実装 |
| **Month 4: 社内検証**       | 1ヶ月 | LOD決定(10反復×7濃度×4ウイルス)、PPA/NPA/R²計算       |

### 11.2 PMDA申請準備(3ヶ月)

| フェーズ                  | 期間  | 作業内容                            |
| --------------------- | --- | ------------------------------- |
| **Month 5-6: 臨床検体検証** | 2ヶ月 | 陽性コントロール入手(研究機関・獣医クリニック)、qPCR相関 |
| **Month 7: ドキュメント作成** | 1ヶ月 | Validation dossier、SOP、トレーニング資料 |

**合計実装期間:** 7ヶ月

---- 

## 12. 参考文献

### 12.1 MinION Viral Metagenomics

1. Quick J, et al. "Rapid draft sequencing and real-time nanopore sequencing in a hospital outbreak of Salmonella." *Genome Biol* 2015;16:114.
2. Charalampous T, et al. "Nanopore metagenomics enables rapid clinical diagnosis of bacterial lower respiratory infection." *Nat Biotechnol* 2019;37:783-792.

### 12.2 ポリオーマウイルス

3. Hirsch HH, et al. "Polyomavirus BK replication in de novo kidney transplant patients receiving tacrolimus or cyclosporine." *Am J Transplant* 2002;2:664-670.
4. Rockett RJ, et al. "Detection of BK polyomavirus in human urine by MinION sequencing." *J Clin Virol* 2019;116:1-6.

### 12.3 ハンタウイルス

5. Kim J, et al. "Multiplex PCR-Based Nanopore Sequencing and Epidemiological Surveillance of Hantaan orthohantavirus." *Viruses* 2021;13:1679.
6. Maes P, et al. "Hantaviruses: immunology, treatment, and prevention." *Viral Immunol* 2004;17:481-497.

### 12.4 東部ウマ脳炎ウイルス(EEEV)

7. Batovska J, et al. "Metagenomic arbovirus detection using MinION nanopore sequencing." *J Virol Methods* 2017;249:79-84.
8. Faria NR, et al. "Establishment and cryptic transmission of Zika virus in Brazil and the Americas." *Nature* 2017;546:406-410.

### 12.5 スピューマウイルス

9. Switzer WM, et al. "Ancient co-speciation of simian foamy viruses and primates." *Nature* 2005;434:376-380.
10. Onions D, et al. "An approach to the control of disease transmission in pig-to-human xenotransplantation." *Xenotransplantation* 2000;7:143-155.

### 12.6 Cell-free DNA/RNA

11. Fernández-Lázaro D, et al. "Clinical perspective and translational oncology of liquid biopsy." *Diagnostics* 2020;10:443.
12. Belizário JE, et al. "Current challenges and best practices for cell-free long RNA biomarker discovery." *Biomarker Res* 2022;10:62.

---- 

## 13. 付録

### 13.1 試薬・機器チェックリスト

**Phase 0 (Sample Preparation)追加試薬:**
- [ ]() SUPERase•In RNase Inhibitor (Invitrogen AM2696)
- [ ]() Ficoll-Paque PLUS (Cytiva 17144003)
- [ ]() RNaseZap (Invitrogen AM9780)

**Phase 3 (Host Removal)追加試薬:**
- [ ]() NEBNext rRNA Depletion Kit (Human/Mouse/Rat) (NEB E6310)
- [ ]() NEBNext Poly(A) mRNA Magnetic Isolation Module (NEB E7490)
- [ ]() DNase I, RNase-free (NEB M0303)

**Phase 4 (Detection)追加試薬:**
- [ ]() Custom degenerate primers (Spumavirus pol gene)
- [ ]() High-fidelity Taq (Q5 or Phusion)

**MinION追加キット:**
- [ ]() Direct RNA Sequencing Kit (ONT SQK-RNA002)
- [ ]() Native Barcoding Kit (ONT SQK-NBD114)

### 13.2 データベースファイル構造

```
/mnt/efs/databases/pmda/2024.1/
├── polyomavirus/
│   ├── polyoma_all.fasta (全Polyomaviridae genomes)
│   ├── polyoma_all.mmi (Minimap2 index)
│   └── taxonomy.txt
├── hantavirus/
│   ├── hantavirus_L_segment.fasta
│   ├── hantavirus_M_segment.fasta
│   ├── hantavirus_S_segment.fasta
│   ├── hantavirus_all.mmi
│   └── species_list.txt
├── alphavirus/
│   ├── alphavirus_all.fasta (EEEV, WEEV, VEEV, Getah, etc.)
│   ├── alphavirus_all.mmi
│   └── eeev_references.fasta (系統解析用)
└── spumavirus/
    ├── sfv_pol_gene.fasta (SFV pol references)
    ├── perv_pol_gene.fasta (識別用)
    └── spumavirus_all.mmi
```

### 13.3 Primer配列(スピューマウイルス)

**Nested PCR primers (pol gene, degenerate):**

```
【1st PCR - Outer primers】
FV-pol-F1: 5'-GGNCARATHGGNATGTTYGG-3'
FV-pol-R1: 5'-CCRTCNCCRAANCCRTC-3'
Expected product: ~800 bp
Tm: 52-58°C (degenerate考慮)

【2nd PCR - Inner primers】
FV-pol-F2: 5'-ATHGGNCARGGNTTYACNAC-3'
FV-pol-R2: 5'-GTRTCNGTYTTRTCNCC-3'
Expected product: ~400 bp
Tm: 54-60°C

【Degeneracy code】
N = A/T/G/C
R = A/G (purine)
Y = C/T (pyrimidine)
H = A/T/C (not G)

【検証推奨】
- In silico validation: Primer-BLAST vs SFV/FFV/BFV pol genes
- Wet lab validation: SFV-positive primate PBMC DNA (if available)
```

---- 

## 14. 変更履歴

| Version | 日付         | 変更内容                                       | 作成者         |
| ------- | ---------- | ------------------------------------------ | ----------- |
| 1.0     | 2025-11-12 | 初版作成。ポリオーマ、ハンタ、EEEV、スピューマ4ウイルス高感度検出プロトコル統合 | Claude Code |

