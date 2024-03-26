# Transcriptome reconstruction for *pseudostellaria heterophylla*
#
# !!IMPORTANT NOTE!!
# ==================
# RNAs in experiment SRX1611055 contain TuMV RNAs.
#
# Steps to reproduce:
# ===================
# 1. install the tools and update the constants under Configurations section
# 2. run this script under target working directory using
#    `rake -f PATH_TO/transcriptome.rake transcriptome:all`
#
require 'fileutils'
require 'rainbow'

namespace :transcriptome do

  # Configurations
  # --------------

  # Path to tools
  SRA_TOOLS_BIN_PATH = '~/sratoolkit-3.0.0/bin'
  FASTQC_BIN_PATH = '~/.miniconda3/envs/toehold/bin'
  TRIMMOMATIC_BIN_PATH = '~/.miniconda3/envs/toehold/bin'
  TRIMMOMATIC_ADAPTERS_PATH = '~/.miniconda3/envs/toehold/share/trimmomatic-0.36-6/adapters'
  SORTMERNA_BIN_PATH = '~/.miniconda3/envs/toehold/bin'
  SORTMERNA_REFS_PATH = '~/data/rRNA-databases'
  TRINITY_BIN_PATH = ENV['TRINITY_HOME']

  # Runs of experiments
  SRXS = {
    # SRX847831: transcriptome of Pseudostellaria CH flower
    # https://www.ncbi.nlm.nih.gov/sra/SRX847831
    'SRX847831' => ['SRR1766410', 'SRR1766414', 'SRR1783710', 'SRR1783711', 'SRR1783712'],

    # SRX1611055: Tanscriptomic analysis of Pseudostellariae Radix using RNA-seq
    # https://www.ncbi.nlm.nih.gov/sra/SRX1611055
    'SRX1611055' => ['SRR3225572']
  }

  # Capacity of the computing environment
  MAX_THREADS = 120
  MAX_MEMORY = '500G'

  # Utility functions
  # -----------------

  @__execute_counter__ = 0
  def execute(one_line, dry_run: false)
    one_line.gsub!(/\s+/, ' ')
    puts format(' %<num>02d  %<cmd>s', num: (@__execute_counter__ += 1), cmd: Rainbow(one_line).yellow)
    system one_line, exception: true unless dry_run
  end

  # Each step as one task
  # ---------------------

  desc 'Download data archives'
  task :prefetch do
    SRXS.values.flatten.each do |id|
      execute "#{SRA_TOOLS_BIN_PATH}/prefetch --progress #{id}"
    end
  end

  desc 'Prepare clean reads'
  task :prepare do
    def fastqc(seqfiles)
      execute "#{FASTQC_BIN_PATH}/fastqc --threads #{MAX_THREADS} #{seqfiles.join(' ')}"
    end

    SRXS.values.flatten.each do |id|
      next if Dir["#{id}/*.clean.fastq"].any?

      # Dump FASTQs from SRAs
      if Dir["#{id}/#{id}*.fastq"].empty?
        execute <<~END_OF_CMD
          #{SRA_TOOLS_BIN_PATH}/fasterq-dump --split-3 --outdir #{id}
            --progress --threads #{MAX_THREADS} #{id}/#{id}.sra"
        END_OF_CMD
        fastqc(Dir["#{id}/#{id}*.fastq"])
      end

      # Trim reads
      if Dir["#{id}/*.trimmed.fastq"].empty?
        if File.exist?("#{id}/#{id}.fastq")
          # Single-end
          execute <<~END_OF_CMD
            #{TRIMMOMATIC_BIN_PATH}/trimmomatic SE -threads #{MAX_THREADS} -phred33
              #{id}/#{id}.fastq #{id}/#{id}.trimmed.fastq
              ILLUMINACLIP:#{TRIMMOMATIC_ADAPTERS_PATH}/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
          END_OF_CMD
        else
          # Pair-end
          execute <<~END_OF_CMD
            #{TRIMMOMATIC_BIN_PATH}/trimmomatic PE -threads #{MAX_THREADS} -phred33
              #{id}/#{id}_1.fastq #{id}/#{id}_2.fastq
              #{id}/#{id}_1.trimmed.fastq #{id}/#{id}_1.trimmed_unpaired.fastq
              #{id}/#{id}_2.trimmed.fastq #{id}/#{id}_2.trimmed_unpaired.fastq
              ILLUMINACLIP:#{TRIMMOMATIC_ADAPTERS_PATH}/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
          END_OF_CMD
        end
        fastqc(Dir["#{id}/*.trimmed.fastq"])
      end

      # Filtering out rRNA
      if Dir["#{id}/*.clean.fastq"].empty?
        if File.exist?("#{id}/#{id}.fastq")
          # Single-end
          execute <<~END_OF_CMD
            #{SORTMERNA_BIN_PATH}/sortmerna #{['-ref'].product(Dir["#{SORTMERNA_REFS_PATH}/*.fasta"]).join(' ')}
              -workdir #{id}/sortmerna/ --threads #{MAX_THREADS} -fastx -other
              -reads #{id}/#{id}.trimmed.fastq
          END_OF_CMD

          FileUtils.mv "#{id}/sortmerna/out/other.fq", "#{id}/#{id}.clean.fastq"
        else
          # Pair-end
          execute <<~END_OF_CMD
            #{SORTMERNA_BIN_PATH}/sortmerna #{['-ref'].product(Dir["#{SORTMERNA_REFS_PATH}/*.fasta"]).join(' ')}
              -workdir #{id}/sortmerna/ --threads #{MAX_THREADS} -fastx -other -paired_out
              -reads #{id}/#{id}_1.trimmed.fastq -reads #{id}/#{id}_2.trimmed.fastq
          END_OF_CMD

          f1 = File.open("#{id}/#{id}_1.clean.fastq", 'w')
          f2 = File.open("#{id}/#{id}_2.clean.fastq", 'w')
          File.open("#{id}/sortmerna/out/other.fq") do |fin|
            until fin.eof?
              4.times { f1 << fin.readline }
              4.times { f2 << fin.readline }
            end
          end
          [f1, f2].each(&:close)
        end
        fastqc(Dir["#{id}/*.clean.fastq"])
      end

      # Remove intermediate files
      Dir["#{id}/*"].each do |file|
        FileUtils.rm_rf(file) unless file =~ /(\.sra|\.clean\.fastq|\.html)$/
      end
    end
  end

  desc 'De novo transcriptome assemble'
  task :assemble do
    SRXS.each do |srx, srrs|
      next if File.exist?("#{srx}_trinity.fasta")

      File.open("#{srx}_trinity_samples", 'w') do |fout|
        srrs.each do |id|
          fout.print("#{srx}\t#{id}\t")
          if File.exist?("#{id}/#{id}.clean.fastq")
            fout.puts "#{id}/#{id}.clean.fastq"
          else
            fout.puts "#{id}/#{id}_1.clean.fastq\t#{id}/#{id}_2.clean.fastq"
          end
        end
      end

      execute <<~END_OF_CMD
        #{TRINITY_BIN_PATH}/Trinity --seqType fq --max_memory #{MAX_MEMORY} --samples_file #{srx}_trinity_samples
          --CPU #{MAX_THREADS} --output #{srx}_trinity --full_cleanup
      END_OF_CMD

      FileUtils.mv "#{srx}_trinity.Trinity.fasta", "#{srx}_trinity.fasta"
    end
  end

  task all: %i[prefetch prepare assemble]
end

task default: 'transcriptome:all'
