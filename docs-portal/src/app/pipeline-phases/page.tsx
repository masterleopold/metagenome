import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { InfoIcon } from "lucide-react";

export default function PipelinePhasesPage() {
  return (
    <div className="container flex-1">
      <div className="flex gap-6">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Core Concepts</Badge>
              <h1 className="text-4xl font-bold mb-4">Pipeline Phases</h1>
              <p className="text-lg text-muted-foreground">
                Six sequential phases for comprehensive pathogen detection and analysis.
              </p>
            </div>

            <Alert className="mb-8">
              <InfoIcon className="h-4 w-4" />
              <div>
                <AlertTitle>Documentation in progress</AlertTitle>
                <AlertDescription>
                  This page is being actively developed. Check back soon for complete documentation.
                </AlertDescription>
              </div>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Overview</h2>
              <p className="text-muted-foreground">
                The MinION analysis pipeline executes six sequential phases, each optimized for
                specific bioinformatics tasks. Each phase runs on dedicated EC2 instances with
                appropriate resources and pre-installed tools.
              </p>
            </section>

        <div className="space-y-6">
          <div className="border-l-4 border-primary pl-6 py-4">
            <div className="flex items-center gap-3 mb-3">
              <span className="bg-primary text-primary-foreground px-3 py-1 rounded-full font-semibold text-sm">
                Phase 1
              </span>
              <h3 className="text-xl font-semibold">Basecalling</h3>
            </div>
            <p className="text-foreground/80 mb-3">
              Converts raw FAST5/POD5 signal data to FASTQ sequences using GPU-accelerated
              Dorado with duplex mode for Q30+ accuracy.
            </p>
            <div className="space-y-2">
              <p className="text-sm font-medium">Key Tools:</p>
              <ul className="list-disc list-inside text-sm text-foreground/80 space-y-1">
                <li>Dorado (GPU basecaller)</li>
                <li>generate_summary_from_fastq.py (post-basecalling metrics)</li>
              </ul>
            </div>
            <div className="mt-3 text-sm text-muted-foreground">
              Instance: g4dn.xlarge (NVIDIA T4 GPU)
            </div>
          </div>

          <div className="border-l-4 border-secondary pl-6 py-4">
            <div className="flex items-center gap-3 mb-3">
              <span className="bg-secondary text-white px-3 py-1 rounded-full font-semibold text-sm">
                Phase 2
              </span>
              <h3 className="text-xl font-semibold">Quality Control</h3>
            </div>
            <p className="text-foreground/80 mb-3">
              Quality assessment using NanoPlot and PycoQC to validate read quality,
              length distribution, and sequencing performance.
            </p>
            <div className="space-y-2">
              <p className="text-sm font-medium">Key Tools:</p>
              <ul className="list-disc list-inside text-sm text-foreground/80 space-y-1">
                <li>NanoPlot (quality visualization)</li>
                <li>PycoQC (comprehensive QC reports)</li>
                <li>qc_check.py (threshold validation)</li>
              </ul>
            </div>
            <div className="mt-3 text-sm text-muted-foreground">
              Instance: c6i.2xlarge
            </div>
          </div>

          <div className="border-l-4 border-amber-500 pl-6 py-4">
            <div className="flex items-center gap-3 mb-3">
              <span className="bg-amber-500 text-white px-3 py-1 rounded-full font-semibold text-sm">
                Phase 3
              </span>
              <h3 className="text-xl font-semibold">Host Removal</h3>
            </div>
            <p className="text-foreground/80 mb-3">
              Aligns reads to Sus scrofa genome using Minimap2 to remove host contamination,
              targeting {'>'}{'>'}90% depletion efficiency.
            </p>
            <div className="space-y-2">
              <p className="text-sm font-medium">Key Tools:</p>
              <ul className="list-disc list-inside text-sm text-foreground/80 space-y-1">
                <li>Minimap2 (long-read aligner)</li>
                <li>Samtools (BAM processing)</li>
                <li>calculate_depletion_rate.py (depletion metrics)</li>
              </ul>
            </div>
            <div className="mt-3 text-sm text-muted-foreground">
              Instance: r6i.2xlarge (memory-optimized)
            </div>
          </div>

          <div className="border-l-4 border-red-500 pl-6 py-4">
            <div className="flex items-center gap-3 mb-3">
              <span className="bg-red-500 text-white px-3 py-1 rounded-full font-semibold text-sm">
                Phase 4
              </span>
              <h3 className="text-xl font-semibold">Pathogen Detection</h3>
            </div>
            <p className="text-foreground/80 mb-3">
              Multi-database pathogen screening using Kraken2, BLAST, and specialized PERV detection
              algorithms for comprehensive coverage of 91 PMDA-designated pathogens.
            </p>
            <div className="space-y-2">
              <p className="text-sm font-medium">Key Tools:</p>
              <ul className="list-disc list-inside text-sm text-foreground/80 space-y-1">
                <li>Kraken2/Bracken (rapid taxonomic classification)</li>
                <li>BLAST (sequence alignment)</li>
                <li>perv_typing.py (PERV-A/B/C subtype identification)</li>
                <li>detect_recombinants.py (PERV recombination analysis)</li>
              </ul>
            </div>
            <div className="mt-3 text-sm text-muted-foreground">
              Instance: c6i.4xlarge
            </div>
          </div>

          <div className="border-l-4 border-green-500 pl-6 py-4">
            <div className="flex items-center gap-3 mb-3">
              <span className="bg-green-500 text-white px-3 py-1 rounded-full font-semibold text-sm">
                Phase 5
              </span>
              <h3 className="text-xl font-semibold">Quantification</h3>
            </div>
            <p className="text-foreground/80 mb-3">
              Absolute quantification of detected pathogens using spike-in normalization
              (PhiX174) to calculate copies per mL plasma.
            </p>
            <div className="space-y-2">
              <p className="text-sm font-medium">Key Tools:</p>
              <ul className="list-disc list-inside text-sm text-foreground/80 space-y-1">
                <li>kraken_quantify.py (Kraken2-based abundance)</li>
                <li>blast_quantify.py (BLAST-based quantification)</li>
                <li>spike_in_normalization.py (PhiX174 normalization)</li>
                <li>absolute_copy_number.py (copies/mL calculation)</li>
              </ul>
            </div>
            <div className="mt-3 text-sm text-muted-foreground">
              Instance: c6i.2xlarge
            </div>
          </div>

          <div className="border-l-4 border-purple-500 pl-6 py-4">
            <div className="flex items-center gap-3 mb-3">
              <span className="bg-purple-500 text-white px-3 py-1 rounded-full font-semibold text-sm">
                Phase 6
              </span>
              <h3 className="text-xl font-semibold">Reporting</h3>
            </div>
            <p className="text-foreground/80 mb-3">
              Generates comprehensive PMDA-compliant reports in PDF, JSON, and HTML formats
              with 91-pathogen checklist and regulatory documentation.
            </p>
            <div className="space-y-2">
              <p className="text-sm font-medium">Key Tools:</p>
              <ul className="list-disc list-inside text-sm text-foreground/80 space-y-1">
                <li>generate_pmda_report.py (multi-format reports)</li>
                <li>generate_pmda_checklist.py (91-pathogen compliance)</li>
                <li>SNS notifications (critical alerts)</li>
              </ul>
            </div>
            <div className="mt-3 text-sm text-muted-foreground">
              Instance: c6i.xlarge
            </div>
          </div>
        </div>

        <section className="space-y-4">
          <h2 className="text-2xl font-semibold">Estimated Timeline</h2>
          <div className="bg-muted p-6 rounded-lg space-y-2">
            <div className="flex justify-between text-sm">
              <span>Phase 1 (Basecalling):</span>
              <span className="font-semibold">2-6 hours</span>
            </div>
            <div className="flex justify-between text-sm">
              <span>Phase 2 (QC):</span>
              <span className="font-semibold">10-20 minutes</span>
            </div>
            <div className="flex justify-between text-sm">
              <span>Phase 3 (Host Removal):</span>
              <span className="font-semibold">30-60 minutes</span>
            </div>
            <div className="flex justify-between text-sm">
              <span>Phase 4 (Pathogen Detection):</span>
              <span className="font-semibold">1-2 hours</span>
            </div>
            <div className="flex justify-between text-sm">
              <span>Phase 5 (Quantification):</span>
              <span className="font-semibold">20-40 minutes</span>
            </div>
            <div className="flex justify-between text-sm">
              <span>Phase 6 (Reporting):</span>
              <span className="font-semibold">5-10 minutes</span>
            </div>
                <div className="border-t pt-2 mt-2 flex justify-between font-semibold">
                  <span>Total Pipeline Duration:</span>
                  <span className="text-primary">4-9 hours</span>
                </div>
              </div>
            </section>
          </div>
        </main>
      </div>
    </div>
  )
}
