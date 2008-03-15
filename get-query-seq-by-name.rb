require 'getoptlong'

opts = GetoptLong.new(
    ["--source", "-s", GetoptLong::REQUIRED_ARGUMENT],
    ["--list", "-l", GetoptLong::OPTIONAL_ARGUMENT],
    ["--name", "-n", GetoptLong::OPTIONAL_ARGUMENT],
 #    ["--qual", "-q", GetoptLong::NO_ARGUMENT],
    ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or !optHash.key?("--source") or (!optHash.key?("--name") and !optHash.key?("--list"))
  $stderr.puts "Usage 1: ruby __.rb -s source.fasta -l query_names.list [-o output.fa ]"
  $stderr.puts "Usage 2: ruby __.rb -s source.fasta -n query_name [-o output.fa ]"
  exit
end

$wanted = {}

if optHash.key?("--list") # list superseed the name argument
  File.new(optHash["--list"], 'r').each do |line|
    if line=~ /^(\S+)/
      $wanted[$1] = 1
    end
  end
elsif optHash.key?("--name")
  $wanted[optHash["--name"]] = 1
end


if optHash.key?("--output")
  outfa = optHash["--output"]
elsif optHash.key?("--list")
  outfa = optHash["--list"] + ".fa"
elsif optHash.key?("--name")
  outfa = optHash["--name"] + ".fa"
end

outfileh = File.new(outfa, "w")

flag = 0
done = {}
File.new(optHash["--source"], "r").each do |line|
  if line=~ /^>(\S+)/
    name = $1
    alter = ''
    if name=~ /^(\S+)\.scf/ # some phrap anormaly 
      alter = $1
    end

    if ( $wanted.key?(name) or $wanted.key?(alter)) and !done.key?(name)
      flag = 1
      done[name] = 1
      outfileh.puts line
    else
      flag = 0
    end
  elsif flag == 1
    outfileh.puts line
  end
end

outfileh.close

squal = optHash["--source"] + ".qual"

if File.exist?(squal) #
  flag = 0
  qualout = File.new(outfa + ".qual", "w")
  
  done = {}
  
  File.new(squal, 'r').each do |line|
    if line=~ /^>(\S+)/
      name = $1
      alter = ''
      if name=~ /^(\S+)\.scf/
        alter = $1
      end
      
      if ( $wanted.key?(name) or $wanted.key?(alter)) and !done.key?(name)
        flag = 1
        done[name] = 1
        qualout.puts line
      else
        flag = 0
      end
    elsif flag == 1
      qualout.puts line
    end
  end
  qualout.close
end

  
