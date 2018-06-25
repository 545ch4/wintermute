require 'json'
require 'json/add/core'
class JSONAble < JSON::GenericObject
  def attributes
    self.class.attributes
  end

  def to_h
    to_hash
  end

  def to_hash
    attributes.inject({}) do |t, key|
      value = send(key)
      value = value.original_sequence if value.kind_of?(Primer)
      t.update key.to_s => value
    end
  end
  
  def self.from_hash(hsh)
    new *(attributes.map{|k| hsh[k.to_sym] || hsh[k.to_s]})
  end
  
  def self.json_create(hsh)
    from_hash(hsh)
  end

  def self.from_json(s)
    from_hash(JSON.parse(s))
  end
  
  def self.json_creatable?
    true
  end

  def to_json(*a)
    to_hash.update('json_class' => self.class.name).to_json(*a)
  end

  def to_s
    self.class.name + ': ' + to_hash.inspect
  end

  def inspect
    to_hash.inspect
  end

end