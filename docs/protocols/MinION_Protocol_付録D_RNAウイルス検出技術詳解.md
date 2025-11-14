# MinION Protocol 付録D: RNAウイルス検出技術詳解

**Version:** 1.0
**作成日:** 2025-11-12
**対象:** ハンタウイルス、東部ウマ脳炎ウイルス(EEEV)等のssRNAウイルス検出

---- 

## 1. RNA分解防止の科学的背景

### 1.1 RNaseの種類と特性

| RNase種類      | 特性                  | 分布                      | 阻害剤                 |
| ------------ | ------------------- | ----------------------- | ------------------- |
| **RNase A**  | エンドヌクレアーゼ、CpNpN↓を切断 | ヒト・ブタ膵臓、血清(\~300 ng/mL) | SUPERase•In, RNasin |
| **RNase H**  | RNA-DNAハイブリッド特異的    | 細胞質                     | EDTA(Mg²⁺要求性)       |
| **RNase 1**  | 非特異的エンドヌクレアーゼ       | 血漿、組織                   | SUPERase•In         |
| **RNase T1** | GpN↓特異的             | 真菌由来、環境汚染               | SUPERase•In         |

### 1.2 血漿中のRNA分解メカニズム

**問題点:**
1. **血液採取時の細胞溶解:** 赤血球・白血球溶解→細胞内RNase放出
2. **血漿中の内在RNase:** RNase 1が\~100 ng/mL存在
3. **環境RNase汚染:** 手指、ピペット、チューブ表面

**半減期データ:**
- 保護なし血漿、室温: t₁/₂ \~15-30分
- EDTA血漿、4°C: t₁/₂ \~2-4時間
- SUPERase•In添加、4°C: t₁/₂ \>24時間
- -80°C凍結: 安定(ただし凍結融解で50%分解)

### 1.3 RNase阻害剤の選択

**SUPERase•In (Invitrogen AM2696) - 推奨**

```
【特性】
- 広スペクトラムRNase阻害剤(RNase A/B/C/1/T1)
- 濃度: 20 U/μL
- 使用濃度: 20 U/mL blood (1 μL/mL)
- 利点: 高活性、RT-PCR阻害なし、-20°C安定

【添加タイミング】
- 血液採取直後(チューブ内に事前添加推奨)
- 10 mL EDTA tube → 200 U (10 μL) pre-aliquot

【コスト】
- 100 U = ¥2,000 → 200 U/sample = ¥4,000
```

**RNasin Plus (Promega N2615) - 代替品**

```
【特性】
- RNase A/B/C阻害
- 濃度: 40 U/μL
- 使用濃度: 1 U/μL in RT reaction
- 利点: RT-PCRでの実績、安価
- 欠点: 血漿中では失活しやすい

【推奨用途】
- RNA抽出後のelute保存
- RT反応時の追加保護
```

### 1.4 RNA保存のベストプラクティス

**短期保存(0-4時間):**
```
EDTA blood → 即座に遠心(4°C) → 血漿を-80°C
- RNase阻害剤必須
- 遠心前の室温放置 <2時間
```

**中期保存(4-24時間):**
```
専用保存チューブ使用
- PAXgene Blood RNA Tube (PreAnalytiX)
- Tempus Blood RNA Tube (Thermo Fisher)
- cf-RNA Preservative Tube (Norgen)

利点: 室温30日間安定
欠点: 高価(¥800-1,200/tube)、DNA抽出に不向き
```

**長期保存(\>1週間):**
```
RNA抽出後、-80°C凍結
- Elution buffer: Nuclease-free water + 1 U/μL RNasin
- Aliquot保存(20 μL × 3本)→ 凍結融解1回限定使用
- 長期(>6ヶ月): RNA later添加(Ambion)
```

---- 

## 2. rRNA除去法の詳細比較

### 2.1 rRNAの問題

**血漿cfRNA組成(健常ヒト、ブタ類似と推定):**
- rRNA: 80-90%
- mRNA: 1-5%
- tRNA: 5-10%
- small RNA (miRNA, etc.): 3-5%
- **ウイルスRNA: \<0.01%**(感染時)

**問題:** MinIONシーケンスの80-90%がrRNA読み→ウイルス検出感度低下

### 2.2 RNase H法(NEBNext rRNA Depletion Kit)

**原理:**
```
1. rRNA特異的DNAプローブ(18S, 28S, 5.8S, 5S)をハイブリダイゼーション
2. RNA-DNAハイブリッド形成
3. RNase HがRNA-DNAハイブリッドのRNA鎖を特異的分解
4. DNase IでDNAプローブ除去
5. rRNA除去RNA回収(mRNA, ウイルスRNA残存)
```

**プロトコル詳細:**

```bash
【Input】Total RNA 1-100 ng

【ステップ1: Probe Hybridization】
- RNA 10 μL + rRNA Depletion Probe (Human/Mouse/Rat) 2 μL
- Hybridization Buffer 6 μL + Nuclease-free Water 2 μL
- Total 20 μL
- 95°C 2分 (RNA変性) → 4°C 30秒 (急冷、プローブアニール)

【ステップ2: RNase H Digestion】
- RNase H (NEB) 5 μL添加
- 37°C 30分インキュベート
- rRNA分解(RNA-DNAハイブリッド特異的)

【ステップ3: DNase I Treatment】
- DNase I 5 μL添加
- 37°C 30分
- DNAプローブ除去

【ステップ4: RNA Cleanup】
- Zymo RNA Clean & Concentrator-5
- Binding Buffer 75 μL + 99.5% EtOH 105 μL添加
- Zymo-Spin ICカラムに添加、12,000×g 30秒
- RNA Wash Buffer 400 μL洗浄×2回
- Nuclease-free Water 15 μL溶出
- → rRNA除去RNA 15 μL (input 10 μLから15 μL回収)

【収量】
- rRNA除去効率: 80-95%
- mRNA/viral RNA回収率: 70-90%
- 最終濃度: input 10 ng → output 1-3 ng (rRNAが大部分)
```

**ブタrRNAへの適用:**

```
【問題】
- NEBNext Human/Mouse/Ratキットはヒト・マウス・ラットrRNA特異的
- ブタ(Sus scrofa) rRNA配列は部分的に相同(哺乳類保存)
- 交差反応率: 推定60-80%(未検証)

【対策】
1. カスタムブタrRNAプローブ設計(IDT xGen Lockdown Probes)
   - ブタ18S rRNA (NR_046261.1)
   - ブタ28S rRNA (NR_046262.1)
   - ブタ5.8S rRNA (NR_046263.1)
   - コスト: ¥300,000-500,000 (初回設計費、100 rxns分)

2. 汎用哺乳類rRNA除去キット
   - QIAGEN FastSelect rRNA Removal Kit
   - Illumina Ribo-Zero Plus (廃番)
```

### 2.3 Magnetic Beads法(Ribo-Zero系)

**原理:**
```
1. ビオチン化rRNA特異的プローブをハイブリダイゼーション
2. Streptavidin磁気ビーズでrRNA-プローブ複合体捕捉
3. 磁気スタンドで分離、上清にmRNA/ウイルスRNA
```

**利点:**
- RNase処理不要(RNA完全長保持)
- 回収率高い(80-95%)

**欠点:**
- キット廃番多い(Illumina Ribo-Zero廃番)
- カスタムプローブ必須(ブタ用)
- コスト高(¥3,000-5,000/sample)

### 2.4 Poly(A)選択との併用

**戦略: 2段階除去**

```
【ステップ1】rRNA除去(RNase H法)
- 80-90% rRNA除去

【ステップ2】Poly(A)選択
- 残存rRNA(poly(A)なし)除去
- poly(A)+ viral RNA濃縮(EEEV等)

【効果】
- rRNA: 98-99%除去(2段階)
- poly(A)+ viral RNA濃縮: 100-1,000倍
- 適用: EEEVに最適、ハンタウイルスはステップ2で失われる
```

---- 

## 3. Poly(A)選択 vs rRNA除去の使い分け

### 3.1 ウイルスゲノムのPoly(A)尾有無

| ウイルス                 | Poly(A)尾 | ゲノム型     | 推奨除去法         |
| -------------------- | -------- | -------- | ------------- |
| **東部ウマ脳炎ウイルス(EEEV)** | **あり**   | ssRNA(+) | **Poly(A)選択** |
| **Getahウイルス**        | あり       | ssRNA(+) | Poly(A)選択     |
| **Chikungunyaウイルス**  | あり       | ssRNA(+) | Poly(A)選択     |
| **ハンタウイルス**          | **なし**   | ssRNA(-) | **rRNA除去**    |
| **インフルエンザウイルス**      | なし       | ssRNA(-) | rRNA除去        |
| **SARSCoV-2**        | あり       | ssRNA(+) | Poly(A)選択     |
| **ノロウイルス**           | あり       | ssRNA(+) | Poly(A)選択     |

### 3.2 Poly(A)選択の詳細プロトコル

**使用キット:** NEBNext Poly(A) mRNA Magnetic Isolation Module (E7490)

```bash
【前処理: DNase処理(必須)】
- Total RNA 50 μL
- DNase I (RNase-free, NEB M0303) 1 U/μg RNA添加
- Protector RNase Inhibitor (Roche) 20 U添加
- 10× DNase I Reaction Buffer 5 μL
- 25°C 30分
- EDTA (50 mM) 5 μL添加
- 65°C 10分(DNase不活化)
- → DNase処理RNA 60 μL

【ステップ1: RNA Denaturation】
- DNase処理RNA 60 μL
- 65°C 5分加熱(2次構造解離、poly(A)尾露出)
- 即座に氷上2分(再アニーリング防止)

【ステップ2: Oligo(dT) Beads Binding】
- 氷上のRNA 60 μLにNEB Oligo(dT) Magnetic Beads 50 μLを添加
- ピペッティング混和10回
- 室温10分インキュベート(ローテーター使用)
- poly(A)+ RNAがoligo(dT) beadsに結合
- rRNA, tRNA, poly(A)- RNAは非結合

【ステップ3: Washing】
- 磁気スタンド(DynaMag-2, Thermo Fisher)で2分静置
- 上清を新しいチューブに保存(option: 2nd round用)
- Beadsに対してWash Buffer (NEBNext kit付属) 200 μL添加
- ピペッティング混和、磁気スタンド2分、上清除去
- Wash 2回繰り返し(計2回洗浄)

【ステップ4: Elution】
- Tris Buffer (10 mM Tris-HCl, pH 7.5) 50 μL添加
- ピペッティング混和
- 65°C 2分加熱(oligo(dT)ハイブリッド解離)
- 即座に磁気スタンドで分離
- 上清(poly(A)+ RNA)を新しいチューブに回収
- → Poly(A)+ RNA 約45 μL

【Optional: 2nd Round Selection】
高純度が必要な場合、上記ステップ2-4を繰り返し
- 1st round eluteを再度oligo(dT) beads処理
- rRNA除去効率: 1 round 95% → 2 rounds 99%

【収量とQC】
- Input total RNA 50 μL (10 ng) → Output poly(A)+ RNA 45 μL (0.5-1 ng)
- rRNA: 98%除去
- poly(A)+ viral RNA回収率: 80-90%
- QC: Agilent RNA 6000 Pico Kit → rRNAピーク消失確認
```

### 3.3 選択フローチャート

```
START: ウイルス種同定
    ↓
    [ウイルスはPoly(A)尾を持つか?]
    ↓                    ↓
   YES                  NO
    ↓                    ↓
Poly(A)選択          rRNA除去
(EEEV, Getah等)     (ハンタウイルス等)
    ↓                    ↓
Direct RNA Kit      cDNA Kit or
or cDNA Kit         Amplicon approach
    ↓                    ↓
MinION Sequencing   MinION Sequencing
```

---- 

## 4. Direct RNA vs cDNA Sequencing比較

### 4.1 技術比較表

| 項目             | Direct RNA (SQK-RNA002)  | Direct cDNA (SQK-DCS109) | PCR-cDNA (LSK114)    |
| -------------- | ------------------------ | ------------------------ | -------------------- |
| **Input量**     | 50 ng poly(A)+ RNA       | 10 ng poly(A)+ RNA       | 1 ng total RNA       |
| **Poly(A)要求**  | 必須(または人工付加)              | 推奨(oligo(dT)プライミング可)     | 不要(random hexamers)  |
| **鎖特異性**       | あり(native strand)        | あり(RT strand)            | なし(dsDNA)            |
| **RNA修飾検出**    | 可能(m⁶A, pseudouridine)   | 不可                       | 不可                   |
| **Read長**      | \~RNA長(最大30 kb報告)        | \~RNA長                   | 増幅bias(短い)           |
| **Throughput** | 低(\~0.5-2 Gb/run)        | 中(\~2-5 Gb/run)          | 高(\~10-20 Gb/run)    |
| **正確度**        | 85-95% (R9.4)            | 90-95%                   | 95-99%               |
| **ウイルス検出感度**   | 低(input量制約)              | 中                        | 高(PCR増幅)             |
| **コスト/run**    | ¥100,000 (kit+flow cell) | ¥100,000                 | ¥80,000 (barcoding込) |

### 4.2 Direct RNA Sequencingのワークフロー

**使用キット:** Oxford Nanopore Direct RNA Sequencing Kit (SQK-RNA002)

**重要注意:** 本キットはpoly(A)尾必須。ハンタウイルスには適用不可(人工poly(A)付加必要)。

```bash
【前提】Poly(A)+ RNA 50 ng以上(濃度 >50 ng/μL推奨)

【ステップ1: Reverse Transcription】
- Poly(A)+ RNA 9 μL (~50 ng)
- RT Adapter (ONT kit付属) 1 μL
- 65°C 5分 → 氷上2分(RT adapter annealing)
- Reverse Transcription Mix:
  - Strand-Switching Primer 1 μL
  - RT Buffer 4 μL
  - RT Enzyme 1 μL
  - Nuclease-free Water 4 μL
- Total 20 μL
- 50°C 10分(cDNA合成)→ 70°C 10分(RT不活化)

【ステップ2: RNA Adapter Ligation】
- RT産物 20 μL
- RNA Adapter (RMX) 2 μL
- Ligation Buffer 10 μL
- T4 DNA Ligase 8 μL
- Total 40 μL
- 室温10分(adapter ligation)

【ステップ3: Cleanup】
- Ligation産物 40 μL
- NEBNext Quick Ligation Reaction Buffer 40 μL添加
- AMPure XP Beads 80 μL添加
- 室温5分ビーズ結合
- 磁気スタンド5分、上清除去
- Wash Buffer (ONT kit) 150 μL洗浄×2回
- Elution Buffer 21 μL溶出
- → Direct RNAライブラリ 21 μL

【ステップ4: Flow Cell Priming】
- 上記セクション6.2参照

【ステップ5: Library Loading】
- Sequencing Buffer (SB) 37.5 μL
- Loading Beads (LB) 25.5 μL
- RNA library 12 μL
- Total 75 μL
- Sample port注入、run開始

【Sequencing条件】
- Chemistry: R9.4.1
- Basecalling model: rna_r9.4.1_70bps_hac.cfg
- Run時間: 24-48時間(低throughputのため長時間推奨)
```

### 4.3 Direct cDNA Sequencing(推奨: EEEV検出)

**使用キット:** Oxford Nanopore Direct cDNA Sequencing Kit (SQK-DCS109)

```bash
【Input】Poly(A)+ RNA 10 ng(Direct RNAより10倍少ない)

【ステップ1: 1st Strand cDNA Synthesis】
- Poly(A)+ RNA 9 μL (10 ng)
- Oligo(dT) VN Primer (ONT kit) 1 μL
- 65°C 5分 → 氷上1分
- Reverse Transcription Mix:
  - Strand-Switching Primer 1 μL
  - RT Buffer 4 μL
  - RT Enzyme Mix 2 μL
  - Nuclease-free Water 3 μL
- Total 20 μL
- 42°C 90分(長時間RT、full-length cDNA合成)
- 70°C 10分(RT不活化)

【ステップ2: 2nd Strand Synthesis】
- 1st strand cDNA 20 μL
- 2nd Strand Synthesis Buffer 10 μL
- 2nd Strand Synthesis Enzyme 5 μL
- Nuclease-free Water 15 μL
- Total 50 μL
- 16°C 10分(2nd strand開始)
- RNase H 1 μL添加
- 37°C 60分(RNA template除去、2nd strand完成)
- 70°C 10分(不活化)

【ステップ3: Cleanup】
- AMPure XP 1.0× (50 μL beads)
- 通常洗浄、25 μL Elution Buffer溶出

【ステップ4: Adapter Ligation】
- dsDNA 24 μL
- Oxford Nanopore Ligation Sequencing Kit (LSK114)使用
- 上記セクション5.1参照

【出力】
- cDNAライブラリ ~15 μL
- 濃度: 50-200 fmol
```

**Direct cDNA利点(EEEVに最適):**
- Direct RNAより10倍少ないinput
- Throughput高い(2-5 Gb/run vs 0.5-2 Gb for Direct RNA)
- Poly(A)+ RNA前処理と組み合わせで高濃縮

---- 

## 5. SMART-9N増幅法の原理と応用

### 5.1 SMART技術の原理

**SMART = Switching Mechanism at 5' End of RNA Template**

```
【メカニズム】
1. Random 9N primerでcDNA合成開始
2. Reverse transcriptaseが5'末端到達
3. Terminal transferase活性で非鋳型C残基3個付加(CCC)
4. SMART oligo (3'-GGG-5')がCCCとハイブリダイゼーション
5. RTがSMART oligoを鋳型として伸長→鋳型switch
6. Full-length cDNAに共通アダプター配列付加(5'/3'両端)

【利点】
- Unbiased全RNA増幅(poly(A)不要)
- 長鎖cDNA合成(最大18.5 kb報告)
- 少量input対応(1 ng - 1 pg)
- ハンタウイルス等poly(A)-ウイルスに有効
```

### 5.2 SMART-9N Protocol

**推奨キット:** Takara SMART-Seq v4 Ultra Low Input RNA Kit

```bash
【Input】Total RNA 10 pg - 10 ng(超高感度)

【ステップ1: 1st Strand cDNA Synthesis】
- RNA 9.5 μL
- 3' SMART-Seq CDS Primer II A (random 9N + adapter) 1 μL
- 72°C 3分(RNA変性)→ 氷上2分
- 5' SMART-Seq v4 Oligonucleotide 1 μL
- First-Strand Buffer 4 μL
- DTT (20 mM) 0.5 μL
- dNTP Mix (10 mM each) 2 μL
- RNase Inhibitor 0.25 μL
- SMARTScribe Reverse Transcriptase 2 μL
- Total 20 μL
- 42°C 90分(cDNA合成 + template switching)
- 70°C 10分(不活化)

【ステップ2: cDNA Amplification】
- 1st strand cDNA 20 μL
- SeqAmp PCR Buffer (2×) 25 μL
- IS PCR Primer (universal primer) 2.5 μL
- Nuclease-free Water 2.5 μL
- SeqAmp DNA Polymerase 1 μL
- Total 50 μL
- PCR: 95°C 1分 → (98°C 10秒 → 65°C 30秒 → 68°C 3分) × 18 cycles → 72°C 10分

【ステップ3: Cleanup】
- AMPure XP 1.0× (50 μL beads)
- 洗浄、40 μL Elution Buffer溶出
- → Amplified cDNA ~38 μL (100-500 ng expected)

【ステップ4: QC】
- Qubit dsDNA HS: 濃度測定
- Agilent High Sensitivity DNA Kit: サイズ分布(500 bp - 10 kb期待)

【ステップ5: MinION Library Prep】
- Amplified cDNA 100-200 fmol
- Oxford Nanopore Ligation Sequencing Kit (LSK114)
- 通常DNAライブラリプロトコル
```

### 5.3 SMART-9Nのハンタウイルス適用

**利点:**
- Poly(A)不要→ハンタウイルス(poly(A)-)検出可能
- L/M/S 3分節同時増幅
- rRNA除去後のlow input (1 ng)から増幅

**欠点:**
- PCR増幅bias(GC-rich領域増幅不良の可能性)
- 18 cyclesでもhost cDNA増幅→ウイルス濃縮効果限定的

**推奨戦略:**
```
rRNA除去 → SMART-9N増幅 → MinION sequencing
or
rRNA除去 → Amplicon RT-PCR(ハンタ特異的) → MinION sequencing (より高感度)
```

---- 

## 6. Amplicon-based検出戦略

### 6.1 ARTIC Network様タイルアンプリコンデザイン

**開発背景:** SARS-CoV-2パンデミックでARTIC Networkが開発した手法

**原理:**
```
1. ウイルスゲノム全体を400 bp amplicons(100 bp overlap)でカバー
2. Primer Pool AとBで交互にamplicons増幅(multiplex PCR干渉回避)
3. Pool A + Pool B → complete genome coverage
4. MinION sequencingでconsensus genome組立
```

### 6.2 ハンタウイルス用タイルアンプリコン設計

**参考論文:** Kim et al., Viruses 2021, "Multiplex PCR-Based Nanopore Sequencing and Epidemiological Surveillance of Hantaan orthohantavirus"

**設計仕様:**

```
【L segment (6,533 nt)】
- 20 amplicons × 400 bp = 8,000 bp coverage(overlap含む)
- Primer Pool A: amplicons 1, 3, 5, 7, 9, 11, 13, 15, 17, 19
- Primer Pool B: amplicons 2, 4, 6, 8, 10, 12, 14, 16, 18, 20

【M segment (3,616 nt)】
- 11 amplicons × 400 bp
- Primer Pool A: 5 amplicons
- Primer Pool B: 6 amplicons

【S segment (1,696 nt)】
- 5 amplicons × 400 bp
- Primer Pool A: 2 amplicons
- Primer Pool B: 3 amplicons

【Total】
- 36 amplicons
- 72 primers (36 forward + 36 reverse)
```

**Primer設計ガイドライン:**

```
1. Conserved region選択
   - Hantaan, Seoul, Dobrava, Puumala virus alignmentから保存領域抽出
   - MAFFT multiple alignment
   - Conservation >80%の領域

2. Primer3でprimer設計
   - Length: 20-25 bp
   - Tm: 60-62°C (±1°C以内)
   - GC%: 40-60%
   - Avoid: secondary structure, self-dimer, hetero-dimer
   - Degenerate base使用可(N, R, Y)最大4 degenerate sites

3. In silico validation
   - Primer-BLAST vs all Hantavirus genomes (NCBI)
   - Off-target check vs Sus scrofa genome

4. Amplicon overlap
   - 100 bp overlap between adjacent amplicons
   - 保証: 1 ampliconがfailしても隣接ampliconでカバー
```

### 6.3 Multiplex RT-PCR Protocol

```bash
【ステップ1: 1st Strand cDNA Synthesis】
- rRNA除去RNA 10 μL
- Random Hexamers (50 μM) 1 μL
- dNTP Mix (10 mM each) 1 μL
- 65°C 5分 → 氷上2分
- 5× RT Buffer 4 μL
- DTT (0.1 M) 1 μL
- RNase Inhibitor (40 U/μL) 1 μL
- SuperScript IV RT (200 U/μL) 1 μL
- Nuclease-free Water 1 μL
- Total 20 μL
- 23°C 10分(random hexamer annealing)
- 50°C 10分(cDNA synthesis)
- 80°C 10分(RT不活化)

【ステップ2: Multiplex PCR (Pool A)】
- cDNA 2.5 μL
- Q5 High-Fidelity 2× Master Mix (NEB M0492) 12.5 μL
- Primer Pool A (各primer 0.015 μM, 18 pairs mixed) 10 μL
- Total 25 μL
- PCR:
  98°C 30秒
  ↓
  (98°C 15秒 → 63°C 5分) × 35 cycles
  ↓
  65°C 5分

【ステップ3: Multiplex PCR (Pool B)】
- 同様にPrimer Pool Bで実施

【ステップ4: Pool Mixing & Cleanup】
- Pool A 25 μL + Pool B 25 μL = 50 μL
- AMPure XP 1.0× (50 μL beads)
- 80% EtOH洗浄×2
- Elution Buffer 20 μL溶出

【ステップ5: QC】
- Agilent High Sensitivity DNA Kit
- Expect: 400 bp broad peak(複数ampliconsの混合)
- Concentration: >10 ng/μL

【ステップ6: Native Barcoding(optional、多検体同時解析)】
- Oxford Nanopore Native Barcoding Kit (SQK-NBD114)
- 各サンプルにbarcode付加(最大24 samples multiplex)

【ステップ7: MinION Sequencing】
- Barcoded amplicons pooling
- Ligation Sequencing Kit (LSK114)でアダプターライゲーション
- 通常シーケンスプロトコル
```

### 6.4 Amplicon vs Metagenomic比較(ハンタウイルス)

| 項目                     | Amplicon-based         | Metagenomic (SMART-9N) |
| ---------------------- | ---------------------- | ---------------------- |
| **感度(LOD)**            | 10-100 copies/mL       | 10⁴-10⁵ copies/mL      |
| **特異度**                | 高(プライマー特異的)            | 中(off-target mapping)  |
| **Genome coverage**    | 高(\>95%, tiled design) | 変動(低力価で部分的)            |
| **Hands-on time**      | 6時間                    | 8時間                    |
| **Cost**               | ¥15,000/sample         | ¥25,000/sample         |
| **Unbiased discovery** | 不可(既知ウイルスのみ)           | 可能(novel variants検出)   |
| **推奨用途**               | 定量・型別、LOD \<100必要時     | 初回スクリーニング、未知ウイルス探索     |

**結論:** ハンタウイルス検出でLOD \<50 copies/mL達成にはAmplicon approach必須。

---- 

## 7. RNA品質評価とトラブルシューティング

### 7.1 RNA Integrity Number (RIN)

**測定機器:** Agilent 2100 Bioanalyzer + RNA 6000 Nano/Pico Kit

**RIN評価:**
```
RIN 10: Intact RNA (18S/28S rRNA峰明瞭、baseline平坦)
RIN 7-9: Partially degraded(使用可能)
RIN 5-6: Moderate degradation(MinION可、qPCR困難)
RIN <5: Severe degradation(使用不可)
```

**血漿cfRNAの特徴:**
- 断片化RNA主体(50-200 nt fragments)
- rRNA混在でRIN測定困難
- rRNA除去後のRIN測定推奨

**目標:**
- rRNA除去前: RIN \>7
- rRNA除去後: RIN測定不可(rRNAピーク消失)、代わりにElectropherogramで評価

### 7.2 QC Checkpoints

| Checkpoint        | 測定項目    | 合格基準                  | 不合格時対策                     |
| ----------------- | ------- | --------------------- | -------------------------- |
| **採血後**           | 処理時間    | \<4時間(4°C)            | cf-RNA Preservative Tube使用 |
| **血漿分離後**         | ヘモリシス   | A414/A540 \<0.2       | 再採血                        |
| **RNA抽出後**        | 濃度      | \>2 ng/μL             | 5 mL血漿使用またはプール             |
|                   | RIN     | \>7                   | RNase混入確認、再抽出              |
| **DNase処理後**      | DNA残存   | qPCR Ct \>35(β-actin) | DNase量2倍、処理時間延長            |
| **rRNA除去後**       | rRNA除去率 | \<10% residual        | 2nd round rRNA除去           |
| **Poly(A)選択後**    | rRNAピーク | 不検出(Bioanalyzer)      | 2nd round poly(A)選択        |
| **Library prep後** | 濃度      | \>50 fmol             | Input増量またはPCR cycles増加     |

### 7.3 トラブルシューティングマトリックス

#### 問題1: RNA収量低い(\<1 ng from 5 mL plasma)

| 原因      | 検証方法           | 対策                              |
| ------- | -------------- | ------------------------------- |
| ヘモリシス   | 血漿の色(赤色)       | 遠心条件改善(1,900×g → 2,500×g)       |
| RNase分解 | RIN \<5        | SUPERase•In 40 U/mL、処理時間短縮      |
| 抽出効率低下  | Carrier RNA未添加 | Carrier RNA (yeast tRNA) 1 μg添加 |
| 血漿量不足   | \<5 mL input   | 10 mL採血                         |

#### 問題2: rRNA除去不十分(\>20% residual)

| 原因          | 検証方法                     | 対策               |
| ----------- | ------------------------ | ---------------- |
| ブタrRNA交差反応低 | Bioanalyzer 18S/28Sピーク残存 | カスタムブタrRNAプローブ設計 |
| RNase H失活   | 陽性コントロール(ヒトRNA)で検証       | 新しいキット使用         |
| RNA 2次構造    | 95°C変性不十分                | 95°C 5分に延長       |

#### 問題3: ウイルスRNA検出されない(qPCR陽性、NGS陰性)

| 原因                 | 検証方法             | 対策                        |
| ------------------ | ---------------- | ------------------------- |
| ウイルス力価低(\<10³)     | qPCR Ct \>35     | Amplicon approach変更       |
| Poly(A)選択でロス       | ウイルスがpoly(A)-    | rRNA除去法に変更                |
| Library prep失敗     | Qubit \<10 fmol  | Input RNA増量、RT条件最適化       |
| Sequencing depth不足 | Total reads \<1M | Run時間延長(48時間)、Flow cell追加 |

#### 問題4: Host RNA混入多(\>90% reads)

| 原因          | 検証方法              | 対策                            |
| ----------- | ----------------- | ----------------------------- |
| rRNA除去不十分   | Bioanalyzer       | 2nd round rRNA除去              |
| Host mRNA混入 | Poly(A)選択でmRNAも濃縮 | Hybridization probe capture追加 |
| ウイルス力価極低    | qPCR Ct \>38      | Amplicon approach必須           |

---- 

## 8. コスト最適化戦略

### 8.1 段階的アプローチ

**Tier 1: 標準スクリーニング(全サンプル)**
```
- Dual DNA/RNA extraction (Zymo Quick-cfDNA/RNA Kit)
- DNA: CpG depletion → MinION (ポリオーマ検出)
- RNA: Poly(A)選択 → cDNA or Direct RNA → MinION (EEEV検出)
- Cost: ¥130,000/sample
```

**Tier 2: 精密検査(Tier 1陽性または疑陽性)**
```
- rRNA除去追加(ハンタウイルス検出)
- Hybridization probe capture(全ウイルス濃縮)
- Cost: +¥15,000 = ¥145,000/sample
```

**Tier 3: 確認試験(Tier 2陽性)**
```
- Amplicon RT-PCR(ハンタウイルス)
- qPCR/qRT-PCR定量
- Sanger確認(スピューマウイルス)
- Cost: +¥5,000 = ¥150,000/sample
```

**結果:** 陽性率5%と仮定
- 95サンプル × ¥130,000 = ¥12,350,000
- 5サンプル × ¥150,000 = ¥750,000
- **Total: ¥13,100,000 / 100サンプル = ¥131,000/sample平均**

従来法¥449,574/sampleと比較し、**71%コスト削減**

---- 

## 9. 参考文献

1. Stark R, et al. "RNA sequencing: the teenage years." *Nat Rev Genet* 2019;20:631-656.
2. Hardwick SA, et al. "Synthetic microbe communities provide internal reference standards for metagenome sequencing and analysis." *Nat Commun* 2018;9:3096.
3. Kim J, et al. "Multiplex PCR-Based Nanopore Sequencing and Epidemiological Surveillance of Hantaan orthohantavirus." *Viruses* 2021;13:1679.
4. Quick J, et al. "Multiplex PCR method for MinION and Illumina sequencing of Zika and other virus genomes directly from clinical samples." *Nat Protoc* 2017;12:1261-1276.
5. Garalde DR, et al. "Highly parallel direct RNA sequencing on an array of nanopores." *Nat Methods* 2018;15:201-206.
6. Belizário JE, et al. "Current challenges and best practices for cell-free long RNA biomarker discovery." *Biomarker Res* 2022;10:62.