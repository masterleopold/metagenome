-- MinION Metagenomics Pipeline Initial Seed Data
-- PMDA 91 Pathogens and Reference Data
-- Created: 2025-01-08

-- ===== PMDA 91 Pathogens Master Data =====
-- Note: This is a subset of the full 91 pathogens list
-- Full list should be imported from official PMDA documentation

BEGIN TRANSACTION;

-- Viral Pathogens (High Priority)
INSERT INTO pmda_pathogens (
    pathogen_code, pathogen_name_ja, pathogen_name_en, pathogen_type,
    risk_level, reference_genome_id, reference_genome_length,
    detection_priority, is_zoonotic, clinical_significance
) VALUES
-- Porcine Endogenous Retroviruses
('PERV-A', 'ブタ内在性レトロウイルスA', 'Porcine endogenous retrovirus A', 'virus',
 'BSL-2', 'NC_003059.1', 8600, 1, true,
 'Critical xenotransplantation safety concern. Potential for human infection.'),

('PERV-B', 'ブタ内在性レトロウイルスB', 'Porcine endogenous retrovirus B', 'virus',
 'BSL-2', 'NC_003060.1', 7900, 1, true,
 'Critical xenotransplantation safety concern. Potential for human infection.'),

('PERV-C', 'ブタ内在性レトロウイルスC', 'Porcine endogenous retrovirus C', 'virus',
 'BSL-2', 'NC_003061.1', 8200, 1, true,
 'Does not infect human cells but may recombine with PERV-A.'),

-- High-risk viral pathogens
('HEV', 'E型肝炎ウイルス', 'Hepatitis E virus', 'virus',
 'BSL-2', 'NC_001434.1', 7176, 1, true,
 'Zoonotic hepatitis. Can cause severe disease in humans.'),

('JEV', '日本脳炎ウイルス', 'Japanese encephalitis virus', 'virus',
 'BSL-3', 'NC_001437.1', 10976, 1, true,
 'Severe neurological disease. Endemic in Asia.'),

('PRRSV', 'ブタ繁殖・呼吸障害症候群ウイルス', 'Porcine reproductive and respiratory syndrome virus', 'virus',
 'BSL-2', 'NC_001961.1', 15428, 2, false,
 'Major economic impact on pig farming. Not zoonotic.'),

('PCV2', 'ブタサーコウイルス2型', 'Porcine circovirus 2', 'virus',
 'BSL-2', 'NC_005148.1', 1768, 2, false,
 'Causes PMWS in pigs. Not known to infect humans.'),

('PRV', '仮性狂犬病ウイルス', 'Pseudorabies virus', 'virus',
 'BSL-2', 'NC_006151.1', 143461, 2, true,
 'Aujeszky disease. Rare human infections reported.'),

('FMDV', '口蹄疫ウイルス', 'Foot-and-mouth disease virus', 'virus',
 'BSL-3', 'NC_011450.1', 8172, 1, false,
 'Highly contagious animal disease. Not directly zoonotic.'),

('ASFV', 'アフリカ豚熱ウイルス', 'African swine fever virus', 'virus',
 'BSL-3', 'NC_044959.2', 189628, 1, false,
 'Severe hemorrhagic fever in pigs. Not zoonotic.'),

('CSFV', '豚熱ウイルス', 'Classical swine fever virus', 'virus',
 'BSL-3', 'NC_002657.1', 12300, 1, false,
 'Highly contagious pig disease. Not zoonotic.'),

('SIV', 'ブタインフルエンザウイルス', 'Swine influenza virus', 'virus',
 'BSL-2', 'NC_026438.1', 13588, 1, true,
 'Can reassort with human influenza. Pandemic potential.'),

('PPV', 'ブタパルボウイルス', 'Porcine parvovirus', 'virus',
 'BSL-2', 'NC_001718.1', 5075, 3, false,
 'Reproductive failure in pigs. Not zoonotic.'),

('EMCV', '脳心筋炎ウイルス', 'Encephalomyocarditis virus', 'virus',
 'BSL-2', 'NC_001479.1', 7835, 3, true,
 'Can cause myocarditis in various species including humans.'),

('RV', '狂犬病ウイルス', 'Rabies virus', 'virus',
 'BSL-3', 'NC_001542.1', 11932, 1, true,
 'Fatal encephalitis. Critical public health concern.'),

('PEDV', 'ブタ流行性下痢ウイルス', 'Porcine epidemic diarrhea virus', 'virus',
 'BSL-2', 'NC_003436.1', 28033, 3, false,
 'Severe diarrhea in pigs. Not zoonotic.'),

('TGEV', 'ブタ伝染性胃腸炎ウイルス', 'Transmissible gastroenteritis virus', 'virus',
 'BSL-2', 'NC_038861.1', 28586, 3, false,
 'Gastrointestinal disease in pigs. Not zoonotic.'),

-- Bacterial Pathogens (High Priority)
('SA', '黄色ブドウ球菌', 'Staphylococcus aureus', 'bacteria',
 'BSL-2', 'NC_007795.1', 2821361, 2, true,
 'Common pathogen. MRSA strains are of particular concern.'),

('SP', '肺炎レンサ球菌', 'Streptococcus pneumoniae', 'bacteria',
 'BSL-2', 'NC_003098.1', 2160842, 2, true,
 'Causes pneumonia, meningitis. Vaccine preventable.'),

('SS', 'ブタレンサ球菌', 'Streptococcus suis', 'bacteria',
 'BSL-2', 'NC_009442.1', 2007491, 1, true,
 'Emerging zoonotic pathogen. Causes meningitis in humans.'),

('EC', '大腸菌', 'Escherichia coli', 'bacteria',
 'BSL-2', 'NC_000913.3', 4641652, 2, true,
 'Pathogenic strains can cause severe disease.'),

('SE', 'サルモネラ', 'Salmonella enterica', 'bacteria',
 'BSL-2', 'NC_003197.2', 4857450, 1, true,
 'Major foodborne pathogen. Multiple serovars.'),

('CT', '破傷風菌', 'Clostridium tetani', 'bacteria',
 'BSL-2', 'NC_004557.1', 2799251, 2, true,
 'Causes tetanus. Vaccine preventable.'),

('CP', 'ウェルシュ菌', 'Clostridium perfringens', 'bacteria',
 'BSL-2', 'NC_008261.1', 3031430, 2, true,
 'Causes gas gangrene and food poisoning.'),

('LA', 'リステリア菌', 'Listeria monocytogenes', 'bacteria',
 'BSL-2', 'NC_003210.1', 2944528, 2, true,
 'Foodborne pathogen. Severe in immunocompromised.'),

('BA', '炭疽菌', 'Bacillus anthracis', 'bacteria',
 'BSL-3', 'NC_007530.2', 5227419, 1, true,
 'Causes anthrax. Bioterrorism concern.'),

('BP', '類鼻疽菌', 'Burkholderia pseudomallei', 'bacteria',
 'BSL-3', 'NC_007434.1', 4074542, 1, true,
 'Causes melioidosis. Endemic in Southeast Asia.'),

('FT', '野兎病菌', 'Francisella tularensis', 'bacteria',
 'BSL-3', 'NC_006570.2', 1892819, 1, true,
 'Causes tularemia. Highly infectious.'),

('YP', 'ペスト菌', 'Yersinia pestis', 'bacteria',
 'BSL-3', 'NC_003143.1', 4653728, 1, true,
 'Causes plague. Historical pandemic pathogen.'),

('MT', '結核菌', 'Mycobacterium tuberculosis', 'bacteria',
 'BSL-3', 'NC_000962.3', 4411532, 1, true,
 'Causes tuberculosis. Major global health concern.'),

('MB', 'ウシ型結核菌', 'Mycobacterium bovis', 'bacteria',
 'BSL-3', 'NC_002945.4', 4345492, 1, true,
 'Causes bovine tuberculosis. Zoonotic.'),

('MA', '鳥型結核菌', 'Mycobacterium avium', 'bacteria',
 'BSL-2', 'NC_008595.1', 5475491, 2, true,
 'Opportunistic pathogen. MAC infections.'),

('CJ', 'カンピロバクター', 'Campylobacter jejuni', 'bacteria',
 'BSL-2', 'NC_002163.1', 1641481, 2, true,
 'Leading cause of bacterial gastroenteritis.'),

('HP', 'ヘリコバクター・ピロリ', 'Helicobacter pylori', 'bacteria',
 'BSL-2', 'NC_000915.1', 1667867, 3, true,
 'Causes peptic ulcers and gastric cancer.'),

('BB', 'ボレリア', 'Borrelia burgdorferi', 'bacteria',
 'BSL-2', 'NC_001318.1', 910724, 2, true,
 'Causes Lyme disease. Tick-borne.'),

('LP', 'レプトスピラ', 'Leptospira interrogans', 'bacteria',
 'BSL-2', 'NC_005823.1', 4627366, 2, true,
 'Causes leptospirosis. Waterborne zoonosis.'),

('TP', '梅毒トレポネーマ', 'Treponema pallidum', 'bacteria',
 'BSL-2', 'NC_000919.1', 1138006, 2, true,
 'Causes syphilis. Sexually transmitted.'),

('CS', 'クラミジア', 'Chlamydia suis', 'bacteria',
 'BSL-2', 'NC_015408.1', 1104821, 3, true,
 'Porcine pathogen. Related to human pathogens.'),

('PM', 'パスツレラ', 'Pasteurella multocida', 'bacteria',
 'BSL-2', 'NC_002663.1', 2257487, 3, true,
 'Causes pasteurellosis. Animal bites.'),

('AP', 'ブタ胸膜肺炎菌', 'Actinobacillus pleuropneumoniae', 'bacteria',
 'BSL-2', 'NC_009053.1', 2274482, 3, false,
 'Causes porcine pleuropneumonia. Not zoonotic.'),

('HPS', 'ヘモフィルス・パラスイス', 'Haemophilus parasuis', 'bacteria',
 'BSL-2', 'NC_011852.1', 2269156, 3, false,
 'Causes Glasser disease in pigs. Not zoonotic.'),

('EP', '豚丹毒菌', 'Erysipelothrix rhusiopathiae', 'bacteria',
 'BSL-2', 'NC_015601.1', 1787941, 2, true,
 'Causes erysipelas in pigs and erysipeloid in humans.'),

('BR', 'ブルセラ', 'Brucella suis', 'bacteria',
 'BSL-3', 'NC_017251.1', 3315175, 1, true,
 'Causes brucellosis. Significant zoonosis.'),

-- Fungal Pathogens
('CN', 'クリプトコッカス', 'Cryptococcus neoformans', 'fungus',
 'BSL-2', 'NC_026745.1', 19050000, 3, true,
 'Opportunistic pathogen. Causes cryptococcosis.'),

('CA', 'カンジダ', 'Candida albicans', 'fungus',
 'BSL-2', 'GCF_000182965.3', 14300000, 3, true,
 'Common opportunistic fungal pathogen.'),

('AF', 'アスペルギルス', 'Aspergillus fumigatus', 'fungus',
 'BSL-2', 'GCF_000002655.1', 29400000, 3, true,
 'Causes aspergillosis in immunocompromised.'),

('PJ', 'ニューモシスチス', 'Pneumocystis jirovecii', 'fungus',
 'BSL-2', 'GCA_001477535.1', 8100000, 3, true,
 'Causes PCP in immunocompromised patients.'),

-- Parasitic Pathogens
('TG', 'トキソプラズマ', 'Toxoplasma gondii', 'parasite',
 'BSL-2', 'NC_001799.1', 80500000, 2, true,
 'Causes toxoplasmosis. Concern for pregnant women.'),

('TC', 'トリヒナ', 'Trichinella spiralis', 'parasite',
 'BSL-2', 'GCA_000181795.1', 64000000, 2, true,
 'Causes trichinellosis. Foodborne parasite.'),

('TS', '有鉤嚢虫', 'Taenia solium', 'parasite',
 'BSL-2', 'GCA_000003025.2', 122000000, 2, true,
 'Causes cysticercosis and taeniasis.'),

('EC-P', 'エキノコックス', 'Echinococcus granulosus', 'parasite',
 'BSL-2', 'GCA_000524195.1', 115000000, 2, true,
 'Causes hydatid disease. Dog-sheep-human cycle.'),

('CP-P', 'クリプトスポリジウム', 'Cryptosporidium parvum', 'parasite',
 'BSL-2', 'NC_006986.1', 9100000, 3, true,
 'Waterborne parasite. Causes cryptosporidiosis.'),

('GD', 'ジアルジア', 'Giardia duodenalis', 'parasite',
 'BSL-2', 'GCA_000002435.1', 11700000, 3, true,
 'Causes giardiasis. Waterborne transmission.'),

('SS-P', '豚回虫', 'Ascaris suum', 'parasite',
 'BSL-1', 'GCA_000298755.1', 273000000, 3, true,
 'Intestinal parasite. Can infect humans.'),

-- Prion Disease
('PRION', 'プリオン病', 'Prion diseases', 'prion',
 'BSL-3', NULL, NULL, 1, true,
 'Transmissible spongiform encephalopathies. No nucleic acid.');

-- ===== Default Analysis Parameters =====
-- These are system-wide default parameters that can be overridden per run

INSERT INTO analysis_parameters (
    run_id, parameter_category, parameter_name, parameter_value
) VALUES
('DEFAULT', 'basecalling', 'model', 'dna_r10.4.1_e8.2_400bps_sup.cfg'),
('DEFAULT', 'basecalling', 'min_qscore', '9'),
('DEFAULT', 'basecalling', 'duplex', 'true'),
('DEFAULT', 'basecalling', 'trim_strategy', 'dna'),

('DEFAULT', 'qc', 'min_read_length', '200'),
('DEFAULT', 'qc', 'max_read_length', '100000'),
('DEFAULT', 'qc', 'min_mean_q', '9'),
('DEFAULT', 'qc', 'min_reads', '100000'),
('DEFAULT', 'qc', 'max_failed_reads_pct', '50'),

('DEFAULT', 'host_removal', 'reference_genome', 'Sus_scrofa_11.1'),
('DEFAULT', 'host_removal', 'min_identity', '90'),
('DEFAULT', 'host_removal', 'min_alignment_length', '100'),

('DEFAULT', 'kraken2', 'confidence_threshold', '0.05'),
('DEFAULT', 'kraken2', 'minimum_hit_groups', '3'),
('DEFAULT', 'kraken2', 'report_zero_counts', 'true'),
('DEFAULT', 'kraken2', 'use_names', 'true'),

('DEFAULT', 'blast', 'evalue_threshold', '1e-10'),
('DEFAULT', 'blast', 'max_target_seqs', '10'),
('DEFAULT', 'blast', 'word_size', '11'),
('DEFAULT', 'blast', 'penalty', '-2'),
('DEFAULT', 'blast', 'reward', '1'),

('DEFAULT', 'assembly', 'min_contig_length', '500'),
('DEFAULT', 'assembly', 'min_coverage', '5'),
('DEFAULT', 'assembly', 'assembler', 'flye'),

('DEFAULT', 'perv', 'min_coverage_breadth', '80'),
('DEFAULT', 'perv', 'min_coverage_depth', '10'),
('DEFAULT', 'perv', 'detect_recombinants', 'true'),
('DEFAULT', 'perv', 'phylogenetic_analysis', 'true'),

('DEFAULT', 'quantification', 'normalization_method', 'spike_in'),
('DEFAULT', 'quantification', 'spike_in_sequence', 'PhiX174'),
('DEFAULT', 'quantification', 'calculate_uncertainty', 'true'),

('DEFAULT', 'reporting', 'format', 'pdf,html,json'),
('DEFAULT', 'reporting', 'include_raw_data', 'false'),
('DEFAULT', 'reporting', 'pmda_checklist', 'true');

-- ===== Test/Development Sample Data =====
-- Only insert if in development environment
-- This should be conditionally executed based on environment

-- Example workflow execution for testing (commented out for production)
/*
INSERT INTO workflow_executions (
    run_id, sample_id, donor_id, status,
    started_at, operator_name, operator_email,
    sample_type, plasma_volume_ml, collection_date,
    flowcell_id, sequencing_kit, pipeline_version
) VALUES
('TEST-RUN-001', 'SAMPLE-001', 'PIG-YUC-001', 'initiated',
 CURRENT_TIMESTAMP, 'Test User', 'test@example.com',
 'plasma', 10.0, '2025-01-08',
 'FAT00001', 'SQK-LSK114', 'v1.0.0');
*/

-- ===== System Configuration =====
-- Store system-wide configuration that doesn't change per run

INSERT INTO analysis_parameters (
    run_id, parameter_category, parameter_name, parameter_value
) VALUES
('SYSTEM', 'retention', 'raw_data_days', '1827'),  -- 5 years for PMDA
('SYSTEM', 'retention', 'analysis_data_days', '1827'),
('SYSTEM', 'retention', 'report_data_days', '3653'),  -- 10 years for reports

('SYSTEM', 'notification', 'critical_pathogen_alert', 'true'),
('SYSTEM', 'notification', 'perv_detection_alert', 'true'),
('SYSTEM', 'notification', 'pipeline_failure_alert', 'true'),
('SYSTEM', 'notification', 'daily_summary', 'true'),

('SYSTEM', 'performance', 'max_parallel_ec2', '5'),
('SYSTEM', 'performance', 'ec2_timeout_minutes', '360'),
('SYSTEM', 'performance', 'lambda_timeout_seconds', '900'),
('SYSTEM', 'performance', 'use_spot_instances', 'false'),

('SYSTEM', 'validation', 'require_positive_control', 'true'),
('SYSTEM', 'validation', 'require_negative_control', 'true'),
('SYSTEM', 'validation', 'min_spike_in_recovery', '50'),
('SYSTEM', 'validation', 'max_spike_in_recovery', '150');

COMMIT;

-- ===== Verification Queries =====
-- Run these to verify the seed data was loaded correctly

-- Count pathogens by type
SELECT pathogen_type, COUNT(*) as count
FROM pmda_pathogens
GROUP BY pathogen_type
ORDER BY count DESC;

-- List high-priority pathogens
SELECT pathogen_code, pathogen_name_en, risk_level, detection_priority
FROM pmda_pathogens
WHERE detection_priority = 1
ORDER BY risk_level DESC, pathogen_name_en;

-- Check system parameters
SELECT parameter_category, COUNT(*) as param_count
FROM analysis_parameters
WHERE run_id IN ('DEFAULT', 'SYSTEM')
GROUP BY parameter_category
ORDER BY parameter_category;