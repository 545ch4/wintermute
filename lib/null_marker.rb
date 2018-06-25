require File.expand_path('json_able.rb', File.dirname(__FILE__))
class NULLMarker < JSONAble
  def initialize
  end

  def call_allele(sequence)
    '?'
  end

  def match(sequence)
    [:match, sequence]
  end

  def self.attributes
    []
  end

  def is_str?
    false
  end
end
