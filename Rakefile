require 'bio'
require 'json'
require 'progressbar'
require 'write_xlsx'

Dir.glob('../lib/tasks/*.rake').each { |r| load r }

class Toehold
  class << self
    attr_accessor :reporter

    def window_search(genome, size: 36, from: 0, to: genome.length, max_padding: 50)
      progressbar = ProgressBar.create(
        title: 'Toehold.window_search',
        total: to - size - from + 1,
        format: '%t: %B | %a | %e'
      )

      from.step(to - size) do |start|
        padding_5p = start > max_padding ? max_padding : start
        padding_3p = (start + size + max_padding) < genome.length ? max_padding : (genome.length - size - start)

        # The context (trigger with padding)
        context = genome[(start - padding_5p)...(start + size + padding_3p)]

        # The toehold
        toehold = new(genome[start...(start + size)])

        yield toehold unless toehold.quantify(context).nil?

        progressbar.increment
      end
    end
  end

  # The toehold trigger RNA from the target genome/sequence that
  # could triggers switches.
  attr_reader :trigger

  # The toehold switch RNA that activates reporter gene expression in
  # response to the corresponding trigger.
  attr_reader :switch, :switch_parts

  def initialize(trigger)
    @trigger = Bio::Sequence::NA.new(trigger).rna

    # Fragment of the switch that complements to the trigger
    sensor = @trigger.reverse_complement

    @switch_parts = [
      # Parts before the sensor
      [
        # A GGG leader sequence that encourages efficient transcription by the T7 RNA polymerase
        case sensor[0...2]
        when /^gg/ then 'g'
        when /^g/ then 'gg'
        else 'ggg'
        end
      ],

      # The sensor
      sensor,

      # Parts between the sensor and the RBS
      [
        'gga',
        'cuuua', # 5-nt stem
        'gaac' # first 4-nt of the loop containing the RBS
      ],

      # The RBS
      'agaggaga',

      # Parts between the RBS and the start codon
      [
        'uaaag' # reverse complement to the 5-nt stem
      ],

      # The start codon
      'aug',

      # Parts after the start condon
      [
        # A 11-nt stem complement to part of sensor
        sensor[25..-1].reverse_complement,
        # One more to shift the coding frame for reporter
        sensor[24],
        # A common 21-nt linker sequence coding for low-molecular-weight amino
        # acids added to the N terminus of the gene of interest
        'aaccuggcggcagcgcaaaag'
      ],

      # The reporter gene
      self.class.reporter
    ].map { |s| Bio::Sequence::NA.new(s.is_a?(Array) ? s.flatten.join : s).rna }

    @switch = Bio::Sequence::NA.new(@switch_parts.join).rna
  end

  def quantify(genome)
    # No in-frame stop codon
    return @quantities = nil if Bio::Sequence::NA.new(switch_parts[6]).translate.include?('*')

    genome = Bio::Sequence::NA.new(genome).rna
    @quantities = JSON.parse(`./quantify.py '#{[genome, trigger, switch].to_json}'`)
  end

  def quantities
    @quantities ||= quantify(trigger)
  end
end

task :coarse_search do
  # Load genomes and genes from FASTA files
  seq = Dir.glob('in/*.fasta').map do |path|
    [File.basename(path, '.fasta'), Bio::FastaFormat.new(File.read(path)).seq]
  end.to_h

  [
    { name: 'BC1', from: 2942, to: 2942 + 509 },
    { name: 'BC2', from: 4390, to: 4390 + 459 }
  ].product(%w[lacZ GFP]).each do |region, reporter|
    File.open("out/#{region[:name]}-#{reporter}.tsv", 'w') do |fout|
      fout.puts [
        'Trigger RNA', 'Switch RNA', 'Favourability',
        'Delta MFE', 'Genome LSS', 'Switch LSS', 'Stemloop NED',
        'Genome MFE', 'Genome MFE Structure',
        'Switch MFE', 'Switch MFE Structure',
        'Complex MFE', 'Complex MFE Structure'
      ].join("\t")

      Toehold.reporter = seq[reporter][0...50]
      Toehold.window_search(seq['TuMV-1'], from: region[:from], to: region[:to]) do |toehold|
        fout.puts [toehold.trigger, toehold.switch, toehold.quantities.values].join("\t")
      end
    end
  end
end

namespace :fine_search do
  task :prerequisites do
    unless Dir['out/*.tsv'].any?
      puts <<~END_OF_DOC
        Coarse search results not found.
        ---
        Run task `coarse_search` first:
        $ rake coarse_search
      END_OF_DOC
      exit
    end

    unless Dir['tmp/*_trinity.fasta'].any?
      puts <<~END_OF_DOC
        Pseudostellaria heterophylla transcriptome not found.
        ---
        Run task `transcriptome:all` first:
        $ rake -f transcriptome.rake transcriptome:all
        Then copy all *_trinity.fasta to
      END_OF_DOC
      exit
    end

    unless File.exist?('tmp/TuMV_genomes.fasta')
      puts <<~END_OF_DOC
        TuMV genomes not found.
        ---
        Download TuMV genomes manually as described in following steps:
        1. search "Turnip mosaic virus" on [NIH Nucleotide database](https://www.ncbi.nlm.nih.gov/nuccore);
        2. select filters including "Species/Viruses" and "Molecule types/genomic DNA/RNA" and "Sequence Type/Nucleotide";
        3. click "send to" and select "complete record" to create file;
        4. move the *sequence.fasta* to *tmp/TuMV_genomes.fasta*.
      END_OF_DOC
      exit
    end

    unless File.exist?('tmp/NC_044183.1.fasta')
      # Download the genome of pseudostellaria heterophylla chloroplast
      # https://www.ncbi.nlm.nih.gov/nuccore/NC_044183.1
      url = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal\&save=file\&db=nuccore\&report=fasta_cds_na\&id=1714607434'
      system("wget -O tmp/NC_044183.1.fasta #{url}")
    end
  end
end

task fine_search: 'fine_search:prerequisites' do
  # Make blast database for TuMV genomes
  system('makeblastdb -in tmp/TuMV_genomes.fasta -dbtype nucl') unless Dir['tmp/TuMV_genomes.fasta.n*'].any?

  # Blast transcriptome to TuMV genomes
  Dir['tmp/*_trinity.fasta'].each do |fasta|
    system <<~END_OF_CMD
      blastn -task megablast -num_threads 4 -db tmp/TuMV_genomes.fasta \
        -query #{fasta} -out #{fasta.gsub('_trinity.fasta', '_blastn_to_TuMV')} \
        -evalue 0.01 -outfmt 6
    END_OF_CMD
  end unless Dir['tmp/*_blastn_to_TuMV'].any?

  # Concatenate all background genomes, except ones mapping to TuMV genomes
  File.open('tmp/backgrounds', 'w') do |fout|
    File.open('tmp/NC_044183.1.fasta').each do |line|
      fout.puts(line =~ /^>/ ? line.gsub('lcl|', '').split[0] : line)
    end

    Dir['tmp/*_trinity.fasta'].each do |fasta|
      exceptions = File.readlines(fasta.gsub('_trinity.fasta', '_blastn_to_TuMV'))
        .map { |line| '>' + line.split("\t").first }.uniq

      File.open(fasta) do |fin|
        line = fin.gets
        loop do
          if line =~ /^>/
            if exceptions.include?(line.split.first)
              nil until (line = fin.gets) =~ /^>/
              next
            end
            line.gsub!('TRINITY', File.basename(fasta, '.fasta'))
          end

          fout.puts line
          break if fin.eof?
          line = fin.gets
        end
      end
    end
  end unless File.exist?('tmp/backgrounds')

  # Make blast database for background genomes
  system('makeblastdb -in tmp/backgrounds -dbtype nucl') unless Dir['tmp/backgrounds.n*'].any?

  # Load coarse search results
  triggers = {}
  %w[BC1 BC2].product(%w[lacZ GFP]).each do |region, reporter|
    triggers[region] ||= {}
    File.open("out/#{region}-#{reporter}.tsv") do |fin|
      fin.readline
      fin.each do |line|
        c = line.chomp.split("\t")
        triggers[region][c[0]] ||= {}
        triggers[region][c[0]]['Switch'] ||= c[1]
        triggers[region][c[0]]['Blast'] ||= 0
        triggers[region][c[0]][reporter] = c[2..-1]
      end
    end
  end

  # Run BLASTN
  triggers.each_key do |region|
    File.open('tmp/query', 'w') do |fout|
      triggers[region].each_key { |trigger| fout << ">#{trigger}\n#{trigger}\n" }
    end

    system('blastn -task blastn -num_threads 4 -db tmp/backgrounds -query tmp/query -out tmp/out -outfmt 6')

    File.open('tmp/out').each do |line|
      c = line.chomp.split("\t")
      triggers[region][c[0]]['Blast'] += Math.log10(c[10].to_f) - 1
    end
  end

  # Output into Excel
  workbook = WriteXLSX.new('out/TuMV.xlsx')
  triggers.each_key do |region|
    worksheet = workbook.add_worksheet(region)
    worksheet.write(0, 0, 'Reporter-independent')
    worksheet.write(0, 3, 'Reporter: lacZ')
    worksheet.write(0, 3 + 2, 'Reporter: GFP')
    worksheet.write(0, 3 + 2 + 2, 'Reporter: lacZ')
    worksheet.write(0, 3 + 2 + 2 + 9, 'Reporter: GFP')
    worksheet.write_row(1, 0, [
      ['Trigger RNA', 'Switch RNA', 'Aggregated Blast E-value'],
      ['Favourability', 'Delta MFE'] * 2,
      [
        'Genome LSS', 'Switch LSS', 'Stemloop NED',
        'Genome MFE', 'Genome MFE Structure',
        'Switch MFE', 'Switch MFE Structure',
        'Complex MFE', 'Complex MFE Structure'
      ] * 2
    ].flatten)

    triggers[region].each_key.with_index do |trigger, index|
      worksheet.write(index + 2, 0, [
        trigger,
        triggers[region][trigger]['Switch'],
        triggers[region][trigger]['Blast'],
        triggers[region][trigger]['lacZ'][0..1],
        triggers[region][trigger]['GFP'][0..1],
        triggers[region][trigger]['lacZ'][2..-1],
        triggers[region][trigger]['GFP'][2..-1]
      ].flatten)
    end
  end
  workbook.close
end

task default: :fine_search
