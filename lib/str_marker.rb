require File.expand_path('json_able.rb', File.dirname(__FILE__))
class STRMarker < JSONAble
  ATTRIBUTES = [:offset, :divider, :chromosome, :core_repeats]
  attr_reader *(ATTRIBUTES)
  attr_writer :stutter_ratios

  def initialize(_offset, _divider, _chromosome, _core_repeats)
    @offset = _offset.to_i
    @divider = _divider.to_i
    @chromosome = _chromosome.to_s
    @core_repeats = _core_repeats
    @stutter_ratios = {}
  end

  def call_allele(sequence)
    ('%i.%i'%([(sequence.length - offset) / divider, (sequence.length - offset) % divider])).sub('.0', '')
  end

  def match(sequence)
    [:match, sequence]
  end

  def self.attributes
    ATTRIBUTES
  end

  def is_str?
    true
  end

  def is_snp?
    false
  end

#  def to_reverse_complement
#    STRMarker.new(@offset, @divider, @chromosome, @core_repeats.map{|core_repeat| SequenceHelper.reverse_complement(core_repeat)})
#  end

  def motif(sequence)
    STRMarker.sequence_to_motif(sequence, core_repeats)
  end

  def self.sequence_to_motif(sequence, _core_repeats = 4)
#    puts "Call sequence_to_motif with _core_repeats => #{_core_repeats.inspect}"

    possible_repeat_sequences = []
    if _core_repeats.kind_of?(Array)
      return sequence if sequence.length <= _core_repeats.first.length
#      puts "There are #{_core_repeats.inspect} core-repeats"
      possible_repeat_sequences = _core_repeats
    else
      repeat_length = _core_repeats
      return sequence if sequence.length <= repeat_length
      possible_repeat_sequences = []
      (sequence.length - repeat_length + 1).times do |start_pos|
        possible_repeat_sequences << sequence[start_pos ... (start_pos + 4)]
      end
      possible_repeat_sequences.uniq!
#      puts "There are #{_core_repeats.inspect} core-repeats"
    end

    repeats = []
    possible_repeat_sequences.each do |possible_repeat_sequence|
      last_i = 0
      while i = sequence[last_i .. -1].index(possible_repeat_sequence)
        repeat = sequence[last_i + i .. -1].match(/(#{possible_repeat_sequence})+/).to_s
        repeats << [last_i + i, repeat.length, last_i + i + repeat.length - 1, repeat, " [#{possible_repeat_sequence}]#{repeat.length / possible_repeat_sequence.length} "] if repeat.length > possible_repeat_sequence.length
        last_i += i + repeat.length
      end
    end
#    puts "Found #{repeats.size} possible repeats"
#    puts "repeats => #{repeats.inspect}"

    motifs = {}
    (1..repeats.size).each do |i|
      combinations = repeats.combination(i).to_a
      combinations.each do |_repeats|
        all_repeats_did_fit = true
        i = 0
        _repeats.each do |_repeat|
          if i <= _repeat[0]
            i = _repeat[2]
          else
            all_repeats_did_fit = false
            break
          end
        end
        if all_repeats_did_fit
          _sequence = sequence.dup
          _repeats.each_with_index do |_repeat, _repeat_index|
            _sequence[_repeat[0] .. _repeat[2]] = (('a'.ord + _repeat_index).chr)*_repeat[1]
          end
          motifs[_repeats] = _sequence
        end
      end
    end
#    puts "motifs => #{motifs.inspect}"

    motif = motifs.sort do |(_repeats_a, _sequence_a), (_repeats_b, _sequence_b)|
      cmp = _sequence_b.chars.inject(0){|t,c| c.ord >= 97 ? t + 1 : t} <=> _sequence_a.chars.inject(0){|t,c| c.ord >= 97 ? t + 1 : t}
      cmp == 0 ? (_repeats_a.size <=> _repeats_b.size) : cmp
    end.first

    return sequence unless motif

    if _core_repeats.kind_of?(Array)
#      puts "motif.first => #{motif.first.inspect}"
#      puts "motif.last (before) => #{motif.last.inspect}"
      motif.last.gsub!(/([A-Z]+)/){|_unmotifed_sequence| STRMarker.sequence_to_motif(_unmotifed_sequence, _core_repeats.first.length)}
#      puts "motif.last (after) => #{motif.last.inspect}"
    else
 #     puts "motif.first => #{motif.first.inspect}"
 #     puts "motif.last => #{motif.last.inspect}"
    end

    _sequence = motif.last.dup
    _repeats = motif.first

    while i = _sequence.index(/[a-z]/)
      _repeat_index = _sequence[i].ord - 97
      _sequence.sub!(/#{_sequence[i]}+/, _repeats[_repeat_index].last)
    end
#    puts "return from sequence_to_motif with _core_repeats => #{_core_repeats.inspect} ==> #{_sequence.strip.gsub(/\s+/, ' ').inspect}"
    return _sequence.strip.gsub(/\s+/, ' ')
  end

  def generate_stutter_motifs(motif, count)
    # TODO @stutter_ratios
    return [] if count < 10
    unmutated_motif = motif.dup.gsub(/\]\d+/) do |motif_part_count|
      if motif_part_count[1..-1].to_i > 2
        "#{motif_part_count}@"
      else
        motif_part_count
      end
    end
    (STRMarker.mutate_motif("#{count} #{unmutated_motif}").flatten.uniq - ["#{count} #{motif}"]).select do |motif_with_count|
      r = false
      if m = motif_with_count.match(/^(\d+) /)
        if m[1].to_i >= 10
          r = true
        end
      end
      r
    end
  end

  def self.mutate_motif(motif_with_count)
    if motif_with_count['@']
      r = [STRMarker.mutate_motif(motif_with_count.dup.sub(/\d+@/){|motif_part_count| motif_part_count[0..-2]})].flatten
      if m = motif_with_count.match(/^(\d+) /)
        if m[1].to_i >= 10
          r += STRMarker.mutate_motif(motif_with_count.dup.sub(/\d+@/){|motif_part_count|  motif_part_count[0..-2].to_i - 2}.sub(/^(\d+)/){|count| (count.to_f * 0.01).round})
          r += STRMarker.mutate_motif(motif_with_count.dup.sub(/\d+@/){|motif_part_count|  motif_part_count[0..-2].to_i - 1}.sub(/^(\d+)/){|count| (count.to_f * 0.1).round})
          r += STRMarker.mutate_motif(motif_with_count.dup.sub(/\d+@/){|motif_part_count|  motif_part_count[0..-2].to_i + 1}.sub(/^(\d+)/){|count| (count.to_f * 0.05).round})
        end
      end
      r.flatten
    else
      [motif_with_count].flatten
    end
  end

end
