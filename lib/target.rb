require File.expand_path('reference_sequences.rb', File.dirname(__FILE__))
require File.expand_path('json_able.rb', File.dirname(__FILE__))
require File.expand_path('primer.rb', File.dirname(__FILE__))
require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))
class Target < JSONAble
  ATTRIBUTES = [:name, :left_flank_reference, :right_flank_reference, :left_flank_iupac, :right_flank_iupac, :orientation, :marker_names, :chromosome]
  attr_reader *(ATTRIBUTES)
  attr_reader :forward_primer, :reverse_primer

  REVERSE = 'r'.freeze
  FORWARD = 'f'.freeze
  MIN_UNIQUE_SEQUENCE_LENGTH = 12
  MAX_UNIQUE_SEQUENCE_LENGTH = 30

  def initialize(name, left_flank_reference, right_flank_reference, left_flank_iupac, right_flank_iupac, orientation, marker_names, chromosome, need_to_find_unique_sequence = true)
    @name = name
    @left_flank_iupac = need_to_find_unique_sequence ? Target.find_shortest_unique_sequence_from_end(left_flank_iupac) : left_flank_iupac
    if need_to_find_unique_sequence && (@left_flank_iupac.nil? || @left_flank_iupac.length > (MAX_UNIQUE_SEQUENCE_LENGTH))
      puts "Something strange happend while shorting left_flank_iupac for #{name}: "
      Target.find_shortest_unique_sequence_from_end(left_flank_iupac, true)
      throw ArgumentError.new("Could not find a shortest unique sequence at the end of '#{left_flank_iupac}' at #{@name}") if @left_flank_iupac.nil?
    end
    @forward_primer = Primer.new(@left_flank_iupac)

    @right_flank_iupac = need_to_find_unique_sequence ? Target.find_shortest_unique_sequence_from_start(right_flank_iupac) : right_flank_iupac
    if need_to_find_unique_sequence && (@right_flank_iupac.nil? || @right_flank_iupac.length > (MAX_UNIQUE_SEQUENCE_LENGTH))
      puts "Something strange happend while shorting right_flank_iupac for #{name}: "
      Target.find_shortest_unique_sequence_from_start(right_flank_iupac, true)
      throw ArgumentError.new("Could not find a shortest unique sequence at the begin of '#{right_flank_iupac}' at #{@name}") if @right_flank_iupac.nil?
    end
    @reverse_primer = Primer.new(@right_flank_iupac)
    
    @orientation = orientation.to_s == REVERSE ? REVERSE : FORWARD
    @marker_names = (marker_names.nil? ? [] : marker_names)
    @chromosome = chromosome
    @enabled = true
  end

  def empty_forward_primer?
    @left_flank_iupac.empty? || @forward_primer.empty?
  end

  def empty_reverse_primer?
    @right_flank_iupac.empty? || @reverse_primer.empty?
  end

  def self.attributes
    ATTRIBUTES
  end

  def unique_id
    @unique_id ||= "#{@left_flank_iupac}|#{right_flank_iupac}"
  end

  def ==(o)
    o.kind_of?(self.class) ? o.unique_id == unique_id : false
  end

  def <=>(o)
    o.kind_of?(self.class) ? o.unique_id <=> unique_id : -1
  end

  def enabled?
    @enabled
  end

  def disabled?
    !@enabled
  end

  def disable!
    @enabled = false
  end

  def enable!
    @enabled = true
  end

  def reverse?
    @orientation == REVERSE
  end

  def forward?
    @orientation == FORWARD
  end

  def to_reverse_complement
    new_reverse = SequenceHelper.reverse_complement(@left_flank_iupac)
    new_forward = SequenceHelper.reverse_complement(@right_flank_iupac)
    Target.new("#{@name}R", new_forward, new_reverse, new_forward, new_reverse, @orientation == REVERSE ? FORWARD : REVERSE, @marker_names.dup, @chromosome, false)
  end

  protected

  def self.find_shortest_unique_sequence(s, start_pos, end_pos, print_debug = false)
    s_as_iupac = SequenceHelper.iupac_to_regex_to_match_iupac(s[(start_pos.kind_of?(Array) ? (-1 * MAX_UNIQUE_SEQUENCE_LENGTH) : 0) .. (end_pos.kind_of?(Array) ? (MAX_UNIQUE_SEQUENCE_LENGTH - 1) : -1)])
    max_count = REFERENCE_SEQUENCES.inject(0){|t, (reference_sequence_name, reference_sequence)| t + reference_sequence.scan(s_as_iupac).size}
#    max_count = 1
    (start_pos.kind_of?(Array) ? start_pos : end_pos).each do |limit|
      s_part = s[(start_pos.kind_of?(Array) ? limit : 0) .. (end_pos.kind_of?(Array) ? limit : -1)]
      part_regexp = SequenceHelper.iupac_to_regex_to_match_iupac(s_part)
      match_count = REFERENCE_SEQUENCES.inject(0){|t, (reference_sequence_name, reference_sequence)| t + reference_sequence.scan(part_regexp).size}
      puts "find_shortest_unique_sequence => #{s_part} => #{part_regexp} => #{match_count}" if print_debug
      return s_part if match_count == max_count && s_part.gsub(/A{2,}/, 'A').gsub(/C{2,}/, 'C').gsub(/G{2,}/, 'G').gsub(/T{2,}/, 'T').size >= MIN_UNIQUE_SEQUENCE_LENGTH
    end
    nil
  end

  def self.find_shortest_unique_sequence_from_start(s, print_debug = false)
    Target.find_shortest_unique_sequence(s, 0, ((MIN_UNIQUE_SEQUENCE_LENGTH - 1) .. s.length).to_a, print_debug)
  end

  def self.find_shortest_unique_sequence_from_end(s, print_debug = false)
    Target.find_shortest_unique_sequence(s, (MIN_UNIQUE_SEQUENCE_LENGTH .. s.length).to_a.map{|x| -1 * x}, -1, print_debug)
  end
end
