require File.expand_path('json_able.rb', File.dirname(__FILE__))
require File.expand_path('primer.rb', File.dirname(__FILE__))
require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))
class Target < JSONAble
  ATTRIBUTES = [:forward_sequence, :reverse_sequence, :orientation, :marker_names]
  attr_reader *(ATTRIBUTES)
  attr_reader :forward_primer, :reverse_primer

  REVERSE = 'r'
  FORWARD = 'f'

  def initialize(_forward_sequence, _reverse_sequence, _orientation, _marker_names = nil)
    @forward_sequence = _forward_sequence.to_s
    @forward_primer = Primer.new(@forward_sequence)
    @reverse_sequence = _reverse_sequence.to_s
    @reverse_primer = Primer.new(SequenceHelper::reverse_complement(@reverse_sequence).to_s)
    @orientation = _orientation.to_s == REVERSE ? REVERSE : FORWARD
    @marker_names = (_marker_names.nil? ? [] : _marker_names)
    @enabled = true
  end

  def empty_forward_primer?
    @forward_sequence.empty? || @forward_primer.empty?
  end

  def empty_reverse_primer?
    @reverse_sequence.empty? || @reverse_primer.empty?
  end

  def self.attributes
    ATTRIBUTES
  end

  def unique_id
    @unique_id ||= "#{@forward_sequence}|#{reverse_sequence}"
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

  def name(i = nil)
    @name ||= (marker_names.join('/') + (i.nil? ? '' : "-#{i}"))
  end
end
