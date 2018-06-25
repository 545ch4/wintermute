#!/usr/bin/env ruby

# Copyright 2018 Sascha Willuweit
#
# Licensed under EUROPEAN UNION PUBLIC LICENCE v. 1.2

require 'rubygems'
ENV['BUNDLE_GEMFILE'] = File.expand_path('Gemfile', File.dirname(__FILE__))
require 'bundler/setup'
require 'pry'
require 'levenshtein'
require 'nokogiri'
require 'fileutils'
require 'json'
require 'parallel'
require 'rubyXL'
require 'optparse'

Dir.glob(File.expand_path('lib/*.rb', File.dirname(__FILE__))).each{|l| require l}

DEBUG = false
VERSION = '0.1'
STATISTICS = FastTargetMatchingAlgorithm::DECISIONS.keys + [:str_detection_success, :str_detection_fail, :str_primers_match, :str_primers_mismatch]

r1_file = ''
sample_name = '<unknown>'
if ARGV.empty?
  ARGV << '-h'
elsif ARGV.last[0] != '-'
  r1_file = ARGV.pop
  unless File.exists?(r1_file)
    puts "R1 file '#{r1_file}' does not exists."
    exit 1
  end
  sample_name = File.basename(r1_file).split('.').first.split('_R1').first
end

options = {:min_reads => 10, :dynamic_q => true, :max_n => 3, :config => "#{File.dirname(__FILE__)}/config/forenseq.json", :min_reads_ratio => 0.01, :min_variant_reads_ratio => 0.05, :statistics => false, :sample_name => sample_name, :no_r2 => false, :primer_trimming => true, :adapter_trimming => true, :calling => true, :survey => false, :match_forward => true, :match_reverse => true, :require_adapter => false, :n_trim => true, :verbose => false, :force => false}
OptionParser.new do |opts|
  opts.banner = "Usage: ./#{File.basename(__FILE__)} [options] <R1 FILE>"

  opts.on("--[no-]calling", "[Do]/[Don't do] STR/SNP calling (default: #{options[:calling].inspect})") do |b|
    options[:calling] = b
  end

  opts.on("--[no-]survey", "[Do]/[Don't] summarize all assigned sequences into one directory/file (default: #{options[:survey].inspect})") do |b|
    options[:survey] = b
  end

  opts.on("--[no-]statistics", "[Do]/[Don't] output a separate statistics file (default: #{options[:statistics].inspect})") do |b|
    options[:statistics] = b
  end

  opts.on("-v", "--[no-]verbose", "Run verbosely (default: #{options[:verbose].inspect})") do |b|
    options[:verbose] = b
  end

  opts.on("-c <filename>", "--config <filename>", String, "Configuarion file for PCR mix used (default: #{options[:config].inspect})") do |filename|
    options[:config] = filename
  end

  opts.on("-f", "--[no-]force", "Overwrite result file(s) (default: #{options[:force].inspect})") do |b|
    options[:force] = b
  end

 opts.on("-o", "--output-calling FILE", String, "STR/SNP calls output filename") do |filename|
    options[:output_calling] = filename
  end

  opts.on("-r", "--references FILE", String, "Assign sequences to references stored in FILE") do |filename|
    options[:references] = filename
  end

  opts.on("--[no-]dynamic-q", "Determine minimal Q-value dynamically based on R1/R2 (default: #{options[:dynamic_q].inspect})") do |b|
    options[:dynamic_q] = b
    options[:min_q] = nil
  end

  opts.on("--[no-]n-trim", "Trim at first N (default: #{options[:n_trim].inspect})") do |b|
    options[:n_trim] = b
  end

  opts.on("--append", "Append results to existing file (requires -o)") do |b|
    options[:append] = b
  end

  opts.on("--no-r2", "Don't automatically determine and load R2 (default: #{options[:no_r2].inspect})") do |b|
    options[:no_r2] = true
  end

  opts.on("--[no-]adapter-trimming", "Do/Don't trim adapter sequences (default: #{options[:adapter_trimming].inspect})") do |b|
    options[:adapter_trimming] = b
  end

  opts.on("--[no-]primer-trimming", "Do/Don't trim primer sequences (default: #{options[:primer_trimming].inspect})") do |b|
    options[:primer_trimming] = b
  end

  opts.on("--require-adapter", "Filter sequences that doesn't meet adapter requirements (specified in config file) (default: #{options[:require_adapter].inspect})") do |b|
    options[:require_adapter] = false
  end

  opts.on("--only x[,y,..,z]", Array, "Only this kind of marker(s): str, snp, x[-chromosomal], y[-chromosomal] and a[utosomal] (default: all)") do |list|
    list.each do |only_item|
      if ['str', 'snp', 'x', 'y', 'a'].include?(only_item.downcase)
        options[("only_#{only_item.downcase}").to_sym] = true
      end
    end
  end

  opts.on("--no-match-forward", "Sequence matching (STR/SNP) is performed only by matchng with reverse primer (default: #{options[:match_forward].inspect})") do
    options[:match_forward] = false
  end

  opts.on("--no-match-reverse", "Sequence matching (STR/SNP) is performed only by matchng with forward primer (default: #{options[:match_reverse].inspect})") do
    options[:match_reverse] = false
  end

  opts.on("--max-n N", Integer, "Maximal number of N within matching primers and R1/R2 first 20 bases (default: #{options[:max_n].inspect})") do |n|
    options[:max_n] = n.to_i
  end

  opts.on("--min-q N", String, "Minimal Q-value for R1/R2 (overrides dynamic)") do |n|
    options[:dynamic_q] = false
    options[:min_q] = n.to_s[0]
  end

  opts.on("--min-reads N", Integer, "Minimal reads to consider a sequence (default: #{options[:min_reads].inspect})") do |n|
    options[:min_reads] = n.to_i
    options[:min_reads_ratio] = nil
  end

  opts.on("--min-reads-ratio N", Float, "Minimal ratio of reads (relative to summarized target reads) to consider a sequence (default: #{options[:min_reads_ratio].inspect})") do |n|
    options[:min_reads_ratio] = n.to_f
    options[:min_reads] = nil
  end

  opts.on("--min-variant-reads N", Integer, "Minimal reads to consider a variant") do |n|
    options[:min_variant_reads] = n.to_i
    options[:min_variant_reads_ratio] = nil
  end

  opts.on("--min-variant-reads-ratio N", Float, "Minimal ratio of reads (relative to summarized target reads with same length) to consider a variant (default: #{options[:min_variant_reads_ratio].inspect})") do |n|
    options[:min_variant_reads_ratio] = n.to_f
    options[:min_variant_reads] = nil
  end

  opts.on("--min-call-ratio N", Float, "Minimal ratio of reads (relative to summarized target reads) to call a sequence") do |n|
    options[:min_call_ratio] = n.to_f
  end

end.parse!

if options[:append] && options[:output_calling].nil?
  puts "No output (-o/--output filename) is given"
  exit 1
end

if options[:calling]
  options[:output_calling] = "#{File.dirname(r1_file)}/#{options[:sample_name]}_CALL.xlsx" unless options[:output_calling]
  if File.exists?(options[:output_calling]) && !options[:append]
    if options[:force]
      FileUtils.rm_f(options[:output_calling])
    else
      puts "Output file '#{options[:output_calling]}' exists already"
      exit 0
    end
  end
end

if options[:survey]
  options[:output_survey] = "#{File.dirname(r1_file)}/#{options[:sample_name]}_SURVEY" unless options[:output_survey]
  if Dir.exists?(options[:output_survey])
    if options[:force]
      FileUtils.rm_rf(options[:output_survey])
    else
      puts "Output directory '#{options[:output_survey]}' exists already"
      exit 0
    end
  end
end

if options[:statistics]
  options[:output_statistics] = "#{File.dirname(r1_file)}/#{options[:sample_name]}_CALL_STATISTICS.xlsx" unless options[:output_statistics]
  if File.exists?(options[:output_statistics])
    if options[:force]
      FileUtils.rm_f(options[:output_statistics])
    else
      puts "Output file '#{options[:output_statistics]}' exists already"
      exit 0
    end
  end
end

r2_file = r1_file.sub('_R1', '_R2')
has_r2_file = !options[:no_r2] && (r1_file != r2_file && File.exists?(r2_file))

puts "WINTERMUTEv#{VERSION} #{DEBUG ? 'in DEBUG mode ' : ''}will process file '#{r1_file}' (with#{'out' unless has_r2_file} R2) using options => #{options.map{|k,v| "#{k}=#{v.inspect}"}.join(', ')}"

print "Reading configuration ... "
unless File.exists?(options[:config])
  puts "Config file '#{options[:config]}' not found."
  exit 1
end
config = JSON.load(File.read(options[:config]))
config['targets'].each do |target|
  if options[:only_str]
    if target.marker_names.inject(true){|mem, marker_name| mem && !config['markers'][marker_name].is_str?}
      target.disable!
    end
  end
  if options[:only_snp]
    if target.marker_names.inject(true){|mem, marker_name| mem && !config['markers'][marker_name].is_snp?}
      target.disable!
    end
  end
  if options[:only_y]
    if target.marker_names.inject(true){|mem, marker_name| mem && config['markers'][marker_name].chromosome != 'Y'}
      target.disable!
    end
  end
  if options[:only_x]
    if target.marker_names.inject(true){|mem, marker_name| mem && config['markers'][marker_name].chromosome != 'X'}
      target.disable!
    end
  end
  if options[:only_a]
    if target.marker_names.inject(true){|mem, marker_name| mem && config['markers'][marker_name].chromosome != 'A'}
      target.disable!
    end
  end
end
puts "Done."
puts "Found #{config['markers'].size} markers defined by #{config['targets'].size} targets => #{config['targets'].select{|t| t.enabled?}.size} targets and #{config['targets'].select{|t| t.enabled?}.map{|t| t.marker_names}.flatten.uniq.size} markers are enabled."

print 'Reading R1 ... '
r1 = FASTQReader.read(r1_file)
puts 'Done.'
r2 = nil
if has_r2_file
  print "Reading R2 ... "
  r2 = FASTQReader.read(r2_file, :reverse_complement => :true)
  puts 'Done.'
  if r1.sequences.size != r2.sequences.size
    puts "R1 and R2 have different sizes (R1:#{r1.sequences.size} vs. R2:#{r2.sequences.size}) -> R2 won't be used"
    r2 = nil
    has_r2_file = false
  end
end

r1_min_q = options[:min_q]
if options[:dynamic_q]
  print 'Determine minimal R1 Q-value ... '
  r1_histogram = [0] * FASTQ::Q_VALUES.length
  r1.qualities.map do |quality| 
    quality[0...10].chars.each do |c|
      r1_histogram[c.ord - FASTQ::MIN_Q_VALUE] += 1
    end
  end
  r1_histogram.each_with_index do |h, i|
    if h > 1
      r1_min_q = (i + FASTQ::MIN_Q_VALUE).chr
      break
    end
#    puts "#{(i + FASTQ::MIN_Q_VALUE).chr} => #{h}"
  end
  puts 'Done.'
end
print "Applying Q-value of #{r1_min_q.inspect} (#{FASTQ.percentify_q(r1_min_q)}%) to R1 ... "
r1.apply_q(r1_min_q)
puts 'Done.'

if has_r2_file
  r2_min_q = options[:min_q]
  if options[:dynamic_q]
    print 'Determine minimal R2 Q-value ... '
    r2_histogram = [0] * FASTQ::Q_VALUES.length
    r2.qualities.map do |quality| 
      quality[0...10].chars.each do |c|
        r2_histogram[c.ord - FASTQ::MIN_Q_VALUE] += 1
      end
    end
    r2_histogram.each_with_index do |h, i|
      if h > 1
        r2_min_q = (i + FASTQ::MIN_Q_VALUE).chr
        break
      end
  #    puts "#{(i + FASTQ::MIN_Q_VALUE).chr} => #{h}"
    end
    puts 'Done.'
  end
  print "Applying Q-value of #{r2_min_q.inspect} (#{FASTQ.percentify_q(r2_min_q)}%) to R2 ... "
  r2.apply_q(r2_min_q)
  puts 'Done.'
end
puts "There are #{r1.sequences.size} sequences to process."

if options[:n_trim]
  print "Trimming at first N ... "
  r1.trim_by_sequence_or_regexp!('N')
  puts "Done. #{r1.length} sequenecs are left."
end

if options[:require_adapter] && config['adapter_trimming'] && config['adapter_trimming']['R1'] && config['adapter_trimming']['R1']
  if config['adapter_trimming']['R1']['type'] == 'by_filename'
    sequence_key = config['adapter_trimming']['R1']['mapping'].keys.select{|key| r1_file.include?(key)}.first
    print "Removing sequences without adapter in R1 (by filename => #{sequence_key}: '#{config['adapter_trimming']['R1']['mapping'][sequence_key]}') ... "
    adapter_regexp = SequenceHelper::string_or_iupac_regexp(config['adapter_trimming']['R1']['mapping'][sequence_key])
    r1.delete_if do |meta, sequence, quality, trim_status|
      sequence.match(adapter_regexp)
    end
    puts "Done. #{r1.length} sequenecs are left."
  end
end

if options[:adapter_trimming] && config['adapter_trimming'] && config['adapter_trimming']['R1'] && config['adapter_trimming']['R1']
  config['adapter_trimming']['R1'].each do |trimming_config|
    if trimming_config['type'] == 'by_filename'
      sequence_key = trimming_config['mapping'].keys.select{|key| r1_file.include?(key)}.first
      print "Trimming adapter in R1 (by filename => #{sequence_key}: '#{trimming_config['mapping'][sequence_key]}') ... "
      count_of_trimmed_sequences = r1.trim_by_sequence_or_regexp!(SequenceHelper::string_or_iupac_regexp(trimming_config['mapping'][sequence_key]))
      puts "Done. Successfully trimmed #{count_of_trimmed_sequences} out of #{r1.sequences.size} sequences (~#{((count_of_trimmed_sequences.to_f * 100.0) / r1.sequences.size.to_f).round}%)."
    elsif trimming_config['type'] == 'general'
      sequence_key = 
      print "Trimming adapter in R1 (general '#{trimming_config['sequence']}') ... "
      count_of_trimmed_sequences = r1.trim_by_sequence_or_regexp!(SequenceHelper::string_or_iupac_regexp(trimming_config['sequence']))
      puts "Done. Successfully trimmed #{count_of_trimmed_sequences} out of #{r1.sequences.size} sequences (~#{((count_of_trimmed_sequences.to_f * 100.0) / r1.sequences.size.to_f).round}%)."
    end
  end
end

joined = []
if has_r2_file
  print 'Joining R1 and R2 ... '
  joined = FASTQ.new(r1, r2)
  puts 'Done.'
else
  joined = r1
end
print "Filtering sequences with more than #{options[:max_n]} Ns at the first 20 bases ... "
joined.delete_if do |meta, sequence, quality, trim_status|
  sequence[0...20].count('N') > options[:max_n]
end
puts "Done. There are #{joined.length} good reads."

print "Assign good reads to targets/markers "
results = Parallel.map(joined.each_slice(1000)) do |joined_slice|
  print '*'

  # prepare hashes
  local_marker_sequences = {}
  config['markers'].keys.each do |marker_name|
    local_marker_sequences[marker_name] = []
  end
  local_target_statistics = {}
  config['targets'].size.times do |i|
    local_target_statistics[i] = STATISTICS.inject({}){|t, x| t.update x => 0}
  end

  algo = FastTargetMatchingAlgorithm.new(1.5, nil, {:max_n => options[:max_n], :scan_for_r2 => has_r2_file, :match_forward => options[:match_forward], :match_reverse => options[:match_reverse], :primer_trimming => options[:primer_trimming]})

  joined_slice.each do |meta, sequence, quality, trim_status|
    config['targets'].each_with_index do |target, i|
      if target.enabled?
        (reason, target_sequence) = algo.match(target, sequence)
        unless target_sequence.nil?
          target.marker_names.each do |marker_name|
            (reason, marker_target_sequence) = config['markers'][marker_name].match(target_sequence)
            unless marker_target_sequence.nil?
              local_marker_sequences[marker_name] << marker_target_sequence
              break
            end
          end
        end
        local_target_statistics[i][reason] += 1
      end
    end
  end
  [local_marker_sequences, local_target_statistics]
end

marker_sequences = {}
config['markers'].keys.each do |marker_name|
  marker_sequences[marker_name] = []
end
target_statistics = {}
config['targets'].size.times do |i|
  target_statistics[i] = STATISTICS.inject({}){|t, x| t.update x => 0}
end

results.each do |local_marker_sequences, local_target_statistics|
  local_marker_sequences.each do |marker_name, sequences|
    marker_sequences[marker_name] += sequences.select{|sequence| !sequence.nil? && sequence.length > 20}
  end
  local_target_statistics.each do |i, local_target_statistic|
    local_target_statistic.each do |reason, count|
      target_statistics[i][reason] += count
    end
  end
end
puts ' Done.'

survey_counts = {}
survey_most_common_sequences = {}
survey_most_common_sequences_tmp = {}
if options[:survey]
  print 'Collecting survey ... '
  survey_counts_bases = ['A', 'C', 'G', 'T', 'N']
  marker_sequences.each do |marker_name, sequences|
    survey_counts[marker_name] = survey_counts_bases.inject({}){|t, base| t.update base => []} unless survey_counts[marker_name]
    sequences.each do |sequence|
      sequence.chars.each_with_index do |base, i|
        if survey_counts[marker_name][survey_counts_bases.first][i].nil?
          survey_counts_bases.each do |_base|
            survey_counts[marker_name][_base][i] = 0
          end
        end
        survey_counts[marker_name][base][i] = survey_counts[marker_name][base][i] + 1 if survey_counts[marker_name][base]
      end
    end
    print '*'
  end
  marker_sequences.each do |marker_name, sequences|
    print '*'
    unless survey_most_common_sequences_tmp[marker_name]
      survey_most_common_sequences_tmp[marker_name] = {}
      sequences.each_with_index do |sequence, i|
        survey_most_common_sequences_tmp[marker_name][sequence[0]] = [] unless survey_most_common_sequences_tmp[marker_name][sequence[0]]
        survey_most_common_sequences_tmp[marker_name][sequence[0]] << i
      end
    end
    could_extend_at_least_one_sequence = true
    next_base_idx = 1
    survey_most_common_sequences[marker_name] = {} unless survey_most_common_sequences[marker_name]
    while could_extend_at_least_one_sequence do
      could_extend_at_least_one_sequence = false
      tmp = {}
      survey_most_common_sequences_tmp[marker_name].each do |s, is|
        if s.length == next_base_idx
          _tmp = {}
          is.each do |i|
            if sequences[i].length > next_base_idx
              next_s = s + sequences[i][next_base_idx]
              _tmp[next_s] = [] unless _tmp[next_s]
              _tmp[next_s] << i
            end
          end
          _tmp_count = _tmp.inject(0){|t, (_tmp_s, _tmp_is)| t + _tmp_is.size}
          if _tmp_count > (is.size * 0.5)
            _tmp.delete_if do |_tmp_s, _tmp_is|
              _tmp_is.size < (sequences.size.to_f * 0.01).to_i
            end
            if _tmp.any?
              could_extend_at_least_one_sequence = true
              tmp.update _tmp
            elsif is.size > (sequences.size.to_f * 0.01).to_i
              survey_most_common_sequences[marker_name][s] = is
            end
          elsif is.size > (sequences.size.to_f * 0.01).to_i
            survey_most_common_sequences[marker_name][s] = is
          end
        end
      end
      survey_most_common_sequences_tmp[marker_name] = tmp
      next_base_idx = next_base_idx + 1
    end
    print '*'
  end

  FileUtils.mkdir_p(options[:output_survey])
  survey_files = {}
  marker_sequences.each do |marker_name, sequences|
    unless survey_files[marker_name]
      survey_files[marker_name] = File.open("#{options[:output_survey]}/#{marker_name}.csv", 'w')
      survey_counts_bases.each do |base|
        survey_files[marker_name].puts "#{base};#{survey_counts[marker_name][base].join(';')}"
      end
      survey_most_common_sequences[marker_name].each do |sequence, is|
        survey_files[marker_name].puts "#{is.size};#{sequence}"
      end
      survey_files[marker_name].puts "count;sequence"
    end
    SequenceHelper::histogram_joiner(sequences).sort{|(a_sequence, a_count), (b_sequence, b_count)| (a_sequence.length == b_sequence.length) ? (b_count <=> a_count) : (a_sequence.length <=> b_sequence.length)}.each do |sequence, count|
      survey_files[marker_name].puts "#{count};#{sequence}"
    end
    print '*'
  end
  survey_files.each do |key, survey_file|
    survey_file.flush
    survey_file.close
  end
  puts ' Done.'
end

if options[:calling]
  unless options[:primer_trimming]
    puts "STR/SNP calling is possible with primer trimming only!"
  else
    references = {}
    if options[:references] && File.exists?(options[:references])
      print 'Reading references ... '
      workbook = RubyXL::Parser.parse(options[:references])
      worksheet = workbook[0]
      worksheet[1..-1].each do |row|
        columns = row.cells.map{|cell| cell.value}
        references[columns[0]] = {} unless references[columns[0]]
        references[columns[0]][columns[1]] = [] unless references[columns[0]][columns[1]]
        references[columns[0]][columns[1]] << columns[6]
      end
      puts 'Done.'
    end
    references_keys = references.keys

    print "Calling alleles and writing sequences to xlsx "
    histogram_joiner_options = options.select{|k, v| [:min_variant_reads, :min_variant_reads_ratio].include?(k)}
    if options[:append] && File.exists?(options[:output_calling])
      workbook = RubyXL::Parser.parse(options[:output_calling])
      worksheet = workbook[0]
      current_row = worksheet.dimension.ref.row_range.last + 1
    else
      workbook = RubyXL::Workbook.new
      worksheet = workbook[0]
      (['Sample', 'Marker', 'Allele', 'Size', 'Reads', 'Ratio'] + (references_keys) + ['Sequence']).each_with_index do |header_cell, i|
        worksheet.add_cell(0, i, header_cell)
      end
      current_row = 1
    end

    marker_sequences.each do |marker_name, sequences|
      alleles = {}
      histogram = SequenceHelper::histogram_joiner(sequences, histogram_joiner_options)
      max_count = histogram.values.max.to_f
      histogram.delete_if do |sequence, count|
        (options[:min_reads] && count < options[:min_reads]) || (options[:min_reads_ratio] && count < (max_count * options[:min_reads_ratio]))
      end
      histogram.sort{|(a_sequence, a_count), (b_sequence, b_count)| (a_sequence.length == b_sequence.length) ? (b_count <=> a_count) : (a_sequence.length <=> b_sequence.length)}.each do |sequence, count|
        allele = config['markers'][marker_name].call_allele(sequence)
        alleles[allele] = [] unless alleles[allele]
        alleles[allele] << [sequence, count]
      end

      alleles.each do |allele, data|
        data.each do |sequence, count|
          ([options[:sample_name], marker_name, allele, sequence.length, count, (count * 100.0) / max_count] + references_keys.map{|reference_sample| (references[reference_sample] && references[reference_sample][marker_name] && references[reference_sample][marker_name].delete(sequence)) ? 'X' : ''} +[sequence]).each_with_index do |cell, i|
            worksheet.add_cell(current_row, i, cell)
          end
          current_row += 1
          print '*'
        end
      end
      references_keys.map do |reference_sample| 
        if references[reference_sample] && references[reference_sample][marker_name]
          references[reference_sample][marker_name].each do |sequence|
            ([options[:sample_name], marker_name, config['markers'][marker_name].call_allele(sequence), sequence.length, 0, 0] + references_keys.map{|reference_sample_x| reference_sample == reference_sample_x ? 'X' : ''} +[sequence]).each_with_index do |cell, i|
              worksheet.add_cell(current_row, i, cell)
            end
            current_row += 1
            print '*'
          end
        end
      end
    end
    workbook.write(options[:output_calling])
    puts ' Done.'
  end
end

if options[:statistics]
  print "Writing statistics to xlsx "
  workbook = RubyXL::Workbook.new
  worksheet = workbook[0]
  ['Sample', 'Target', 'Reason', 'Count'].each_with_index do |header_cell, i|
    worksheet.add_cell(0, i, header_cell)
  end
  current_row = 1
  target_statistics.each do |i, statistics|
    STATISTICS.each do |reason|
      if config['targets'][i].enabled?
        [options[:sample_name], config['targets'][i].name(i), reason.to_s, statistics[reason]].each_with_index do |cell, i|
          worksheet.add_cell(current_row, i, cell)
        end
        current_row += 1
      end
    end
  end
  workbook.write(options[:output_statistics])
  puts ' Done.'
end
