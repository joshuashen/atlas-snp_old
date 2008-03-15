require 'getoptlong'

opts = GetoptLong.new(
    ["--source", "-s", GetoptLong::OPTIONAL_ARGUMENT],
    ["--length", "-l", GetoptLong::REQUIRED_ARGUMENT],
    ["--prefix", "-p", GetoptLong::OPTIONAL_ARGUMENT],
    ["--qual", "-q", GetoptLong::NO_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or !optHash.key?("--length")
  $stderr.puts "Usage: ruby __.rb -s source.fasta -l length_of_each_batch [-p prefix_of_output] [-q]"
  exit
end

if !optHash["--prefix"]
  prefix = optHash["--source"] + "_batch"
else
  prefix = optHash["--prefix"]
end

$ll = optHash["--length"].to_i

landmarks = {}
count = 0
num = 1
out = prefix + "_#{num}.fa" 

outh = File.new(out, "w")

if !optHash["--source"]
  inputf = $stdin
else
  inputf = File.new(optHash["--source"], "r")
end

inputf.each do |line|
  if line=~ /^>(\S+)/
    name = $1
    if count > $ll
      landmarks[name] = 1
      outh.close
      count = 0
      num += 1
      out = prefix + "_#{num}.fa"
      outh = File.new(out, "w")
    end
    
    outh.puts line
  else
    count += line.chomp.length
    outh.puts line
  end
end
outh.close

inputf.close

if !optHash["--qual"] or !optHash["--source"]
  exit
end

## 

squal = optHash["--source"] + ".qual"

if !File.exist?(squal)
  if (optHash["--source"] =~ /^(\S+)\.fa$/ or optHash["--source"] =~ /^(\S+)\.fna$/)
    squal = $1 + ".qual"
  end
end

exit unless File.exist?(squal)

num = 1
out = prefix + "_#{num}.fa.qual"
outh = File.new(out, "w")
File.new(squal, 'r').each do |line|
  if line=~ /^>(\S+)/
    name = $1
    if landmarks.key?(name)
      outh.close
      
      num += 1
      out = prefix + "_#{num}.fa.qual"
      outh = File.new(out, "w")
    end

    outh.puts line
  else
    outh.puts line
  end
end

outh.close

      
