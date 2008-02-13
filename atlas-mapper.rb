# do blat and cross_match

require 'getoptlong'

opts = GetoptLong.new(
    ["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT],
    ["--query", "-q", GetoptLong::REQUIRED_ARGUMENT],
#    ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
    ["--minscore", "-m", GetoptLong::OPTIONAL_ARGUMENT],
		["--cutoff", "-c", GetoptLong::OPTIONAL_ARGUMENT],
		["--blat", "-b", GetoptLong::OPTIONAL_ARGUMENT],
		["--crossmatch", "-x", GetoptLong::OPTIONAL_ARGUMENT],
    ["--oneOff", "-n", GetoptLong::OPTIONAL_ARGUMENT],
    ["--minIdentity", "-i", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT],
    ["--xmOnly", "-z", GetoptLong::NO_ARGUMENT],
    ["--short", "-s", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or !optHash.key?("--reference") or !optHash.key?("--query")
  $stderr.puts "Usage: ruby __.rb -q query.fasta -r reference.fasta [-m minscore -c cutoff] [-b path_to_blat -x path_to_cross_match]  [-n 1 -i minIdentity_blat] [-z] [--short]"
  $stderr.puts "  note: cutoff is the percentage for best hit cutoff, default 0.995"
  exit
end

fa = optHash["--query"]
reffa = optHash["--reference"]

if optHash.key?("--oocCut")
	oocCut = optHash["--oocCut"].to_i
else
	oocCut = 1024
end

if optHash.key?("--blat")
	blat = File.expand_path(optHash["--blat"])
else
	blat = "blat"
end

if optHash.key?("--crossmatch")
	crossMatch = File.expand_path(optHash["--crossmatch"])
else
	crossMatch = "cross_match"
end

if optHash.key?("--minscore")
	$minscore = optHash["--minscore"].to_i
else
	$minscore = 24
end

if optHash.key?("--cutoff")
	$cutoff = optHash["--cutoff"].to_f
else
  $cutoff = 0.995
end 


if optHash.key?("--oneOff")
  $oneOff = optHash["--oneOff"].to_i
else
  $oneOff = 0
end

if optHash.key?("--minIdentity")
  $minIdentity = optHash["--minIdentity"].to_f
else
  $minIdentity = 90
end


kmer = 11
$unitSize = 100000

$index = {}
$reflist = {}
$outlist  = {}

tempDir = File.dirname(File.expand_path(fa)) + "/Mapping_of_" + File.basename(fa)
package = File.dirname(File.expand_path($0))

# get absolute path for the reference file
absref = File.expand_path(reffa) 
absquery = File.expand_path(fa)
xmoutfile = File.expand_path(tempDir) + '/cross_match.vs.ref'


envDir = absref + '.Env4mapping/'

# read meta information from ref env
refdivs = []
File.new(envDir+"meta.info").each do |line|
  if line=~ /^pieceLength\s+(\d+)/
    $unitSize = $1.to_i
  elsif line=~ /^division\s+(\S+)/
    refdivs << $1
  end
end

system("mkdir #{tempDir}") unless File.directory?(tempDir)

Dir.chdir(tempDir)
blatPsl = "blat.psl"
oocfile = envDir + "11mer.ooc"
bestHits = blatPsl + ".best_hits"

if !optHash.key?("--xmOnly")
  # do BLAT

  $stderr.puts "\nDoing BLAT"

  
  # options: 1. ooc, 2. oneOff, 3. minIdentity
  ## -oneOff=0 default, 

  
  # clean dir
  system("rm -rf #{blatPsl}") if File.exist?(blatPsl)
  system("rm -rf #{blatPsl}.gz") if File.exist?(blatPsl + ".gz")
  

  
  # refs.each do |refDivFa|
  refdivs.each do |refdiv|
    refDivFa = envDir + "ref_divisions/" + refdiv
    #   refDivFa.chomp!
    pslf = "temp.psl"
    cmd = "#{blat} #{refDivFa} #{absquery} -ooc=#{oocfile} -oneOff=#{$oneOff} -minIdentity=#{$minIdentity} #{pslf}"
    
    $stderr.puts cmd
    system(cmd)
    system("cat #{pslf} >> #{blatPsl} ")
    system("rm -rf #{pslf}")
  end

  # get best hit
  $stderr.puts "\nPicking best hits"
  cmd = "ruby #{package}/pick_best_hits_from_BLAT_psl.rb #{blatPsl} #{$cutoff} > #{bestHits}"
  system(cmd)
  
# could be forked to a new process
  if !File.exist?(blatPsl + '.gz')
    Process.fork { system("gzip #{blatPsl}"); }
  end
end


# reading the ref environment 
refindexFile = envDir + 'fragments.index'
File.new(refindexFile, 'r').each do |line|
  cols = line.split(/\s+/)
  refName,s,e,layer,ff = cols[0],cols[1].to_i,cols[2].to_i,cols[3],cols[4]
  $index[refName] = {} unless $index.key?(refName)
  # puts "#{refName}\t#{s}\t#{e}\t#{layer}\t#{ff}"
  $index[refName][s] = layer + '/' + ff
end


# # group reads 
# according to BLAT psl output
$readMap = {}
File.new(bestHits,"r").each do |line|
  cols = line.split(/\s+/)
  next if cols.size < 11 
  flag = cols[-1]
  next if flag == 'repeat'
  readName, refName, tstart = cols[0], cols[2], cols[7].to_i
  next unless $index.key?(refName)
  
  start = (tstart / $unitSize ).to_i * $unitSize + 1
  
  reffile = envDir + 'ref_pieces/' + $index[refName][start]
  
  readfile = File.basename(reffile) + '.reads'
  $reflist[readfile] = reffile
  $readMap[readName] = readfile
  $outlist[readfile] = ''  unless $outlist.key?(readfile)

end

xmDir = "Temp4xm" 

system("rm -rf #{xmDir}") if File.exist?(xmDir)

system("mkdir #{xmDir}") 
xmDirabs = File.expand_path(xmDir)
Dir.chdir(xmDir)

# do cross_match
$stderr.puts "\nDoing cross_match"
xmo = File.new(xmoutfile, 'w')

def readFasta(ff, str)
  flag =0
  readfile = ''
  
  # distribute fasta files
  File.new(ff, 'r').each do |line|
    if line=~ /^>(\S+)/
      name = $1
      if $readMap.key?(name)
        readfile = $readMap[name] 
        flag = 1
        str[readfile] << line
      else
        flag = 0
      end
    elsif flag == 1
      str[readfile] << line
    end
  end
end

readFasta(absquery, $outlist)

$outlist.each_key do |readfile|
  faname = readfile+'.fa'
#   $stderr.puts faname
  faout = File.new(faname, 'w')
  faout.puts $outlist[readfile]
  faout.close
  $outlist[readfile] = ''
end

# distribution quality files

qualfile = absquery + '.qual'
## $stderr.puts qualfile
if File.exist?(qualfile)
  readFasta(qualfile, $outlist)
  
  $outlist.each_key do |readfile|
    faout = File.new(readfile+'.fa.qual', 'w')
    faout.puts $outlist[readfile]
    faout.close
    $outlist[readfile] = ''
  end
end

# str = ''
$outlist.each_key do |readfile|

  readf4shell = readfile.gsub('|', '\|')
  
  reffile = $reflist[readfile].gsub('|', '\|')

  if optHash.key?("--short")
    # for short reads
    #  -bandwidth 6 -gap_init -2 -penalty -1 -gap_ext -1 -raw  -masklevel 101
    xm = `#{crossMatch} #{readf4shell}.fa #{reffile} -minscore #{$minscore} -bandwidth 6 -gap_init -2 -penalty -1 -gap_ext -1 -raw  -masklevel 101 -discrep_lists  2>/dev/null`
  else
    xm = `#{crossMatch} #{readf4shell}.fa #{reffile} -minscore #{$minscore} -raw -discrep_lists  2>/dev/null`

  end


  xmo.puts xm.scan(/\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\(\d+\).*|^\s[S|D|I]\S*\s+\d+\s+\S+\(\d+\).*/)

#  str << xm.scan(/\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\(\d+\).*|^\s[S|D|I]\S*\s+\d+\s+\S+\(\d+\).*/).to_s + "\n"

end
#xmo.puts str
xmo.close
system("rm -rf #{xmDirabs}")
