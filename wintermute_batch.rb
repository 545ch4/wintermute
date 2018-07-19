#!/usr/bin/env ruby

# Copyright 2018 Sascha Willuweit
#
# Licensed under EUROPEAN UNION PUBLIC LICENCE v. 1.2

files = []
new_argv = ARGV.reverse
new_argv.delete_if do |file|
	break if file.start_with?('-')
	if file.include?('fastq') && !file.include?('.json')
		files << file if File.exists?(file)
		true
	else
		false
	end
end
new_argv.reverse!

puts "Processing #{files.size} FASTQ files"
files.each do |file|
	cmd_line = "#{"#{__FILE__}".sub('_batch', '')} #{new_argv.join(' ')} #{file}"
	system(cmd_line)
end
