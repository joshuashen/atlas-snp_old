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
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or !optHash.key?("--reference") or !optHash.key?("--query")
  $stderr.puts "Usage: ruby __.rb -q query.fasta -r reference.fasta [-m minscore -c cutoff] [-b path_to_blat -x path_to_cross_match]  [-n 1 -i minIdentity_blat]"
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
system("mkdir #{tempDir}") unless File.directory?(tempDir)

Dir.chdir(tempDir)

# do BLAT
blatPsl = "blat.psl"
$stderr.puts "\nDoing BLAT"
oocfile = envDir + "11mer.ooc"

# options: 1. ooc, 2. oneOff, 3. minIdentity
## -oneOff=0 default, 

refs = `ls #{envDir}ref_divisions/div*fa`


# clean dir
system("rm -rf #{blatPsl}") if File.exist?(blatPsl)
system("rm -rf #{blatPsl}.gz") if File.exist?(blatPsl + ".gz")

if refs.length < 1 
  refs = absref
end

refs.each do |refDivFa|
  refDivFa.chomp!
  pslf = "temp.psl"
  cmd = "#{blat} #{refDivFa} #{absquery} -ooc=#{oocfile} -oneOff=#{$oneOff} -minIdentity=#{$minIdentity} #{pslf}"

  $stderr.puts cmd
  system(cmd)
  system("cat #{pslf} >> #{blatPsl} ")
  system("rm -rf #{pslf}")
end

# get best hit
bestHits = blatPsl + ".best_hits"
$stderr.puts "\nPicking best hits"
cmd = "ruby #{package}/pick_best_hits_from_BLAT_psl.rb #{blatPsl} #{$cutoff} > #{bestHits}"
system(cmd)

# could be forked to a new process
if !File.exist?(blatPsl + '.gz')
  Process.fork { system("gzip #{blatPsl}"); }
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
  
  if $outlist.key?(readfile)
    $outlist[readfile] << "\n#{readName}"
  else
    $outlist[readfile] = readName
  end
end

xmDir = "Temp4xm" 
system("mkdir #{xmDir}") unless File.exist?(xmDir)
xmDirabs = File.expand_path(xmDir)
Dir.chdir(xmDir)

# do cross_match
$stderr.puts "\nDoing cross_match"
xmo = File.new(xmoutfile, 'w')


$outlist.each_key do |readfile|
  oo = File.new(readfile, 'w')
  oo.puts $outlist[readfile]
  oo.close
  # get read fasta and quality
  readf4shell = readfile.gsub('|', '\|')
  
  cmd = "perl #{package}/get-query-seq-by-name.pl -s #{absquery} -l #{readf4shell} -q"
  system(cmd)
  # do cross_match
  reffile = $reflist[readfile].gsub('|', '\|')
  
  #	$stderr.puts " #{readf4shell} vs #{reffile}"
  xm = `#{crossMatch} #{readf4shell}.fa #{reffile} -minscore #{$minscore} -raw -discrep_lists  2>/dev/null`
#	$stderr.puts "#{crossMatch} #{readf4shell}.fa #{reffile} -minscore #{$minscore} -raw -discrep_lists  2>/dev/null"
  xmo.puts xm.scan(/\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\(\d+\).*|^\s[S|D|I]\S*\s+\d+\s+\S+\(\d+\).*/)
end

xmo.close
system("rm -rf #{xmDirabs}")
