module SequenceHelper
  SEQUENCE_REGEXP = {
    'A' => '[NA]',
    'G' => '[NG]',
    'C' => '[NC]',
    'T' => '[NT]',
    'N' => 'N'
  }

  IUPAC_REGEXP = {
    'A' => '[NA]',
    'G' => '[NG]',
    'C' => '[NC]',
    'T' => '[NT]',
    'Y' => '[NCT]',
    'R' => '[NAG]',
    'S' => '[NGC]',
    'W' => '[NAT]',
    'K' => '[NTG]',
    'M' => '[NAC]',
    'B' => '[NSKYCGT]',
    'D' => '[NRKWAGT]',
    'H' => '[NMYWACT]',
    'V' => '[NMSRACG]',
    'N' => '.',
    'a' => 'A?',
    'c' => 'C?',
    'g' => 'G?',
    't' => 'T?',
    'y' => '[CT]?',
    'r' => '[AG]?',
    's' => '[GC]?',
    'w' => '[AT]?',
    'k' => '[TG]?',
    'm' => '[AC]?',
    'b' => '[CGT]?',
    'd' => '[AGT]?',
    'h' => '[ACT]?',
    'v' => '[ACG]?',
    'n' => '[ACGT]?'
  }
  
  COMPLEMENT_MAP = {
   'A' => 'T',
   'T' => 'A',
   'U' => 'A',
   'G' => 'C',
   'C' => 'G',
   'Y' => 'R',
   'R' => 'Y',
   'S' => 'S',
   'W' => 'W',
   'K' => 'M',
   'M' => 'K',
   'B' => 'V',
   'D' => 'H',
   'H' => 'D',
   'V' => 'B',
   'N' => 'N'
  }

  def self.complement(s)
   s.to_s.chars.map{|c| COMPLEMENT_MAP[c] || '.'}.join
  end

  def self.reverse_complement(s)
   complement(s).reverse
  end

  def self.histogram(arr)
    hsh = {}
    arr.each do |elem|
      elem
      hsh[elem] = hsh[elem].to_i + 1
    end
    hsh
  end

  def self.n_to_nucleotide_map(sequences, _options = {})
    options = {:min_variant_reads => 1}.update(_options)
    map = []

    nucleotide_counts = sequences.first.chars.inject([]){|t, x| t << {}}
    sequences.each do |elem|
      elem.chars.each_with_index do |c, i|
        nucleotide_counts[i][c] = nucleotide_counts[i][c].to_i + 1
      end
    end

    nucleotide_counts.each_with_index do |hsh, i|
      max_count = hsh.values.max.to_f
      hsh.delete_if do |c, count|
        c != 'N' && ((options[:min_variant_reads] && count < options[:min_variant_reads]) || (options[:min_variant_reads_ratio] && count < (max_count * options[:min_variant_reads_ratio])))
      end
      if hsh.size == 2 && hsh.keys.include?('N')
        map << [i, (hsh.keys - ['N']).first]
      end
    end

    map
  end

  def self.histogram_joiner(arr, _options = {})
    options = {:min_variant_reads => 1}.update(_options)
    arr_by_length = {}
    arr.each do |elem|
      arr_by_length[elem.length] = [] unless arr_by_length[elem.length]
      arr_by_length[elem.length] << elem
    end

    hsh = {}
    arr_by_length.each do |l, elems|
      n_to_nucleotide_map(elems, options).each do |i, c|
        elems.size.times do |j|
          elems[j][i] = c if elems[j][i] == 'N'
        end
      end

      elems.each do |elem|
        hsh[elem] = hsh[elem].to_i + 1
      end
    end
    hsh
  end

  def self.sequence_to_regexp(s)
    Regexp.new(s.chars.map{|c| SEQUENCE_REGEXP.has_key?(c) ? SEQUENCE_REGEXP[c] : c}.join)
  end

  def self.string_or_iupac_regexp(s)
    s.nil? || s.length == 0 ? '' : (IUPAC_REGEXP.keys.inject(false){|t, iupac_char| t || s.include?(iupac_char)} ? iupac_to_regex(s) : s)
  end

  def self.iupac_to_regex(s)
    Regexp.new(s.chars.map{|c| IUPAC_REGEXP.has_key?(c) ? IUPAC_REGEXP[c] : c}.join)
  end

  def self.iupac_to_rregex(s)
    Regexp.new(s.chars.map{|c| IUPAC_REGEXP.has_key?(c) ? IUPAC_REGEXP[c] : c}.join + '.*$')
  end
end