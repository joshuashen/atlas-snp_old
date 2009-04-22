require 'getoptlong'

#Program Overview:
#     'atlasSubstitution.rb' is used to parse out the substitutions from the crossmatch result and list
#  all of these substitutions including their corresponding coordination on reference, coverage, number of
#  variation, the variant reads information etc. As a necessary step, this program plays the important
#  role of collecting necessary information from mapping result and doing preparation for final SNP calling.
#     This program takes four required parameters, first parameter is the crossmatch file introduced by
#  an '--crossmatch' option, the second option '--reference' is associate with the fasta file of the reference genome.
#  The third parameter is an option '--nqs', followed by the '_.NQS' file. The option '--output' requires
#  users to denote the path of the output file including a string as the prefix of the output file.
#  After running this program, the result will be saved in the file named '[prefix].SNP.list'. 
#    This program also provide several other optional parameters, which are tunable. For example,
#  it will filter out those substitutions found in such a read, whose 'maximum substitutions'
#  and 'maximum indels' is greater than 5.
#    
#Command Line Arguments:    
#  Usage: ruby atlasSubstitution.rb -x [cross_match] -r [reference] -n [_.NQS] -o [outputPrefixPath] [-s max_substitution -g max_gap -a adjust_quality_slope -c]
#  -x is a required option, which must be appended by the full path of crossmatch file.
#  -r is a required option. User must denote the full path of the reference genome.
#  -n is a required option. The full path of the _.NQS file must be appended.
#  -o is a required option. User must provide the path including the prefix of the output file.
#  -a is an optional argument which specifies the linear adjustment slope for quality scores. Default 0.13; 0 means no adjustment.
#  -s is an optional argument which sets the maximum number of substitutions contained in a variant read to be a true variant read.
#  -g is an optional argument which sets the maximum number of indels contained in a variant read to be a true variant read.
#  -c is an optional argument. If it is appended, the program will calculate the coverage depth of mapping.
#
#Example:
#  1. Call substitution from the crossmatch file 'chr2.xm.part.0' using reference 'chr2' and 'reads.part_0.NQS' as the NQS file.
#  Print the substitution list into file 'part_0.xm.SNP.list'.The maximum number of subsitutions and indels a read can have to
#  be a variant read are set as default value, 5.
#      ruby atlasSubstitution_V3.0.rb -x chr2.xm.part.0 -r chr2 -o part_0.xm -n reads.part_0.NQS
#  2. Call substitution from the crossmatch file 'chr2.xm.part.0' using reference 'chr2' and 'reads.part_0.NQS' as the NQS file.
#  Print the substitution list into file 'part_0.xm.SNP.list'. The maximum number of subsitutions and indels a read can have to
#  be a variant read are set to 3.
#      ruby atlasSubstitution_op_V3.1.rb -x chr2.xm.part.0 -r chr2 -o part_0.xm -n reads.part_0.NQS -s 3 -g 3
class AtlasSubs

#Initialize options
def AtlasSubs.processArguments()
opts = GetoptLong.new(
    ["--crossmatch", "-x", GetoptLong::REQUIRED_ARGUMENT],
    ["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT],
    ["--output", "-o", GetoptLong::REQUIRED_ARGUMENT],
    ["--nqs","-n", GetoptLong::REQUIRED_ARGUMENT],
    ["--coverage", "-c", GetoptLong::OPTIONAL_ARGUMENT],
    ["--maxsub", "-s", GetoptLong::OPTIONAL_ARGUMENT],
    ["--maxindel", "-g", GetoptLong::OPTIONAL_ARGUMENT],
    ["--slope", "-a", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

AtlasSubs.usage() if (optHash.key?("--help"));

AtlasSubs.usage() if (optHash.empty?);
return optHash

end

#Usage information
def AtlasSubs.usage(msg='')
    unless (msg.empty?)
        puts "\n#{msg}\n"
    end
    puts "
Program Overview:
     'atlasSubstitution.rb' is used to parse out the substitutions from the crossmatch result and list
  all of these substitutions including their corresponding coordination on reference, coverage, number of
  variation, the variant reads information etc. As a necessary step, this program plays the important
  role of collecting necessary information from mapping result and doing preparation for final SNP calling.
     This program takes four required parameters, first parameter is the crossmatch file introduced by
  an '--crossmatch' option, the second option '--reference' is associate with the fasta file of the reference genome.
  The third parameter is an option '--nqs', followed by the '_.NQS' file. The option '--output' requires
  users to denote the path of the output file including a string as the prefix of the output file.
  After running this program, the result will be saved in the file named '[prefix].SNP.list'. 
    This program also provide several other optional parameters, which are tunable. For example,
  it will filter out those substitutions found in such a read, whose 'maximum substitutions'
  and 'maximum indels' is greater than 5.
    
Command Line Arguments:    
  Usage: ruby atlasSubstitution.rb -x [cross_match] -r [reference] -n [_.NQS] -o [outputPrefixPath] [-s max_substitution -g max_gap -a adjust_quality_slope -c]
  -x is a required option, which must be appended by the full path of crossmatch file.
  -r is a required option. User must denote the full path of the reference genome.
  -n is a required option. The full path of the _.NQS file must be appended.
  -o is a required option. User must provide the path including the prefix of the output file.
  -a is an optional argument which specifies the linear adjustment slope for quality scores. Default 0.13; 0 means no adjustment.
  -s is an optional argument which sets the maximum number of substitutions contained in a variant read to be a true variant read.
  -g is an optional argument which sets the maximum number of indels contained in a variant read to be a true variant read.
  -c is an optional argument. If it is appended, the program will calculate the coverage depth of mapping.

Example:
  1. Call substitution from the crossmatch file 'chr2.xm.part.0' using reference 'chr2' and 'reads.part_0.NQS' as the NQS file.
  Print the substitution list into file 'part_0.xm.SNP.list'.The maximum number of subsitutions and indels a read can have to
  be a variant read are set as default value, 5.
      ruby atlasSubstitution_V3.0.rb -x chr2.xm.part.0 -r chr2 -o part_0.xm -n reads.part_0.NQS
      
  2. Call substitution from the crossmatch file 'chr2.xm.part.0' using reference 'chr2' and 'reads.part_0.NQS' as the NQS file.
  Print the substitution list into file 'part_0.xm.SNP.list'. The maximum number of subsitutions and indels a read can have to
  be a variant read are set to 3.
      ruby atlasSubstitution_op_V3.1.rb -x chr2.xm.part.0 -r chr2 -o part_0.xm -n reads.part_0.NQS -s 3 -g 3"
  exit(2);
end

#Initialize the object of class AtlasSubstitution
def initialize(optHash)
    @optHash=optHash
    setParameters()
end

#Initialize the parameters
def setParameters()

if (!@optHash.key?('--crossmatch') or !@optHash.key?("--reference") or !@optHash.key?("--output") or !@optHash.key?("--nqs")) then
        AtlasSubs.usage("Option missing!")
        exit(2);
    end

if @optHash.key?("--maxsub")
  @maxsub = @optHash["--maxsub"].to_f
else
  @maxsub = 5.0
end

if @optHash.key?("--maxindel")
  @maxindel = @optHash["--maxindel"].to_f
else	
  @maxindel = 5.0
end

if @optHash.key?("--slope")
  @slope = @optHash["--slope"].to_f
else
  @slope = 0.13
end

if @optHash.key?("--coverage")
    @coverage=true
else
    @coverage=false
end

@outputFileName=File.basename(@optHash["--output"])
@outputFilePath=File.dirname(@optHash["--output"])
@outputFile=@optHash["--output"]
@nqsFile=@optHash["--nqs"]
@xmFile=@optHash["--crossmatch"]
@refFile=@optHash["--reference"]

end

$basechange={"A"=>"T","T"=>"A","G"=>"C","C"=>"G","N"=>"N","X"=>"X","*"=>"*"}

# cross_match result pattern
$pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/
###Creat hash for the NQS file#####

#Hash .NQS file
def hashNQS()
    @hash_nqs={}
    nqsReader=File.open(@nqsFile,"r")
    nqsReader.each do |line|
        cols=line.split(/\s+/)
    
        name=cols[0]
        next if name=="readName"
        length=cols[1]
        dist=cols[2]
        qual=cols[3].to_i
        pass=cols[4]
    
        str_result=length+'.'+pass
        str=name+'.'+dist
        if @hash_nqs[str]==nil
            @hash_nqs[str]=str_result
        end
    end
    nqsReader.close
    $stderr.puts @hash_nqs.size
end

#Compute coverage, determine "Swap", "SNP" or "MNP"
def compute(name, ref, span, snps)
  return if span.length < 1

  span.sort! {|a,b| a[0] <=> b[0]}
  head = span.shift
  ss,ee = head[0],head[1]
  array = []
  
  array << ss 
  span.each do |breaks|
    array << breaks[0]
    array << breaks[1]
  end
  array << ee 

  while array.size > 0
    s = array.shift
    e = array.shift
#    $stderr.puts "#{ref}\t#{s}\t#{e}"
    (s..e).each do |i|
      @coverage[ref][i] += 1 #Compute coverage on base i
    end
  end

  snps.each_key do |pos|
    refbase = @seq[ref][pos-1,1].upcase
    curbase = snps[pos][:snpbase]
    if snps.key?(pos + 1) or snps.key?(pos + 2) or snps.key?(pos - 1) or snps.key?(pos - 2)
      if snps.key?(pos + 1) and refbase == snps[pos+1][:snpbase] and curbase == @seq[ref][pos,1].upcase
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos - 1) and refbase == snps[pos-1][:snpbase] and curbase == @seq[ref][pos-2,1].upcase 
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos + 2) and refbase == snps[pos+2][:snpbase] and curbase == @seq[ref][pos+1,1].upcase
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos - 2) and refbase == snps[pos-2][:snpbase] and curbase == @seq[ref][pos-3,1].upcase
        snps[pos][:info] << "swap;"
      elsif snps.key?(pos + 1) or snps.key?(pos - 1) 
        snps[pos][:info] << "mnp;"
      else
        snps[pos][:info] << "snp;"        
      end
    else
      snps[pos][:info] << "snp;"      
    end
    @snp[ref][pos] = '' unless @snp[ref].key?(pos)
    @snp[ref][pos] << snps[pos][:info]
  end
end

# count the occurence of homo-polymer runs in a short sequence
def homocount(str)
  s = str.size
  lasts = ''
  homo = {}
  st = 0
  flag = 0
  0.upto(s-1) do |i|
    cur = str[i,1].upcase
    if cur == lasts
      if flag == 0 # this is the second 
        homo[i-1] = 2
        st = i - 1
      else # this is the third or more
        homo[st] += 1
      end
      flag = 1
    else
      flag = 0
    end
    lasts = cur
  end
  if homo.size > 0
    top = homo.keys.sort {|a,b| homo[b] <=> homo[a]}[0]
    xx = homo[top]
  else
    xx = 1
  end
  return xx
end

#Initialize reference genome
def initializeReference
    ref = ''
    @seq={}
    @snp={}
    refReader=File.open(@refFile,"r")
    refReader.each do |line|
        if line=~ /^>(\S+)/
          ref = $1
          @seq[ref] = ''
          @snp[ref] = {}
        else
          @seq[ref] << line.chomp
        end
    end
    refReader.close
    $stderr.puts @seq.size
    $stderr.puts @snp.size
end

#Initialize coverage
def initializeCoverage()
    @coverage={}
    @seq.each_key do |ref|
      @coverage[ref] = Array.new(@seq[ref].size + 1) {|i|  0}
    end
end

#Parse substitutions
def parseSubs()
    offset = 0
    @query,@ref,qsize,score,dir = '','', 0, 0,''
    flag = 0
    @span = []
    @snps = {}
    sub, gap, tail = 0, 0, 0
    c1=0
    $stderr.puts "start to read xm file!"
    xmReader=File.open(@xmFile,"r")
#Read crossmatch file
    xmReader.each do |line|
        if line.match($pattern)
            
          compute(@query,@ref,@span,@snps)
          score,sub,del,ins,@query,qstart,qend,qright,tstrand,target,d1,d2,d3,lab =$1.to_f,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11,$12,$13,$14

          @span  = []
          @snps = {}
          @ref = ''
            qsize = qend + qright
            dir = '+'
            ratio = (qend - qstart + 1)/qsize.to_f
            gap = del + ins
            tail = [qstart, qright].max
            if sub > @maxsub or gap > @maxindel 
              flag = 0
              next
            else
              flag = 1
            end
            if target=~ /^(\S+)\_(\d+)\_(\d+)$/
              @ref = $1
            offset = $2.to_i - 1
            elsif target=~ /^(\S+)$/
              @ref = $1
            offset = 0
            end
    
            if d1 =~ /\((\d+)\)/
              tstart, tend = d2.to_i, d3.to_i
            elsif d3=~ /\((\d+)\)/
              tstart, tend  = d1.to_i, d2.to_i 
            end
            if tstart > tend
              tstart, tend = tend, tstart  # wicked 
              dir = '-'
            end
            s = tstart + offset
            e = tend + offset

    #    compute(query,s,e)
            @span << [s,e]
    
        elsif flag > 0 and line=~ /^\s(\S+)\s+(\d+)\s+(\S+)\((\d+)\)\s+(\d+)\s+(\S+)/
          type, qplace, tplace = $1, $2.to_i, $5.to_i
          ii = $3
          qual = $4
          env = $6
    
          if dir == '-' 
            ii = $basechange[ii]
          end
          if type =~ /^D/
      # deletion on read      
             block = 1
            if type =~ /^D\-(\d+)/
                block = $1.to_i
            end
      
            tstart = tplace
            if dir == '-'  #
               tstart = tstart - block + 1
             end
      
            tend = tstart + block
            @span << [tstart + offset, tend + offset]

          elsif type =~ /^S/ # substitution
            tstart = tplace + offset
            next if ii == 'N'
            dist = qsize - qplace 
            @snps[tstart] = {}
            @snps[tstart][:snpbase] = ii
            str=@query+'.'+dist.to_s
      
            if @hash_nqs[str]!=nil
                result=@hash_nqs[str]
                if result=~ /(\S+)\.(\S+)/
                    length=$1
                    pass=$2
                end
            else
                $stderr.puts @query
                $stderr.puts tplace
                $stderr.puts dist
                $stderr.puts "####"
                c1+=1
            end
            @snps[tstart][:info] = "#{ii}(#{qual})#{@query}(#{dist})(#{score}/#{qsize})#{dir}#{env}(#{sub}/#{gap}/#{tail})(#{pass}/#{length})"
          end
        else
            $stderr.puts line
        end
    end
    xmReader.close
    compute(@query,@ref,@span, @snps)
    
    $stderr.puts c1
end

#Print result
def printSubs()
 @refEnv=''
 snpWriter = File.new(@outputFile+".SNP.list", 'w')
#snpout.puts "refName\tcoordinate\trefBase\thomopolymer\trefEnv\tcoverage\tSNPBase\tadjustedQual\toriQual\tnumSNPReads\tnumAlterReads\treads_info"
 snpWriter.print "refName\tcoordinate\trefBase\thomopolymer\trefEnv\tcoverage\tSNPBase\tadjustedQual\toriQual\tnumVariantReads\tnumAlternativeReads\treads_info\n"
 @snp.keys.sort.each do |ref|
  @snp[ref].keys.sort.each do |pos|
    bases = {}
    t = 0
    @snp[ref][pos].split(';').each do |r|
      if r=~ /^(\S)\((\d+)\)(\S+)\((\d+)\)\((\S+)\/(\d+)\)/
        base = $1
        rname = $3
        t += 1
        if !bases.key?(base)
          bases[base] = {}
          bases[base][:adjQual] = 0
          bases[base][:oriQual] = 0
          bases[base][:num] = 0
        end
        bases[base][:num] += 1
        phredQual = $2.to_i
        bases[base][:oriQual] += phredQual
        dist = $4.to_i
        qsize = $6.to_i
        
        if dist < 100 
          adjqual1 = [phredQual - ( (100 - dist) * @slope ),0].max
        else
          adjqual1 = phredQual
        end
        bases[base][:adjQual] += adjqual1
      end
    end
    array = bases.keys.sort {|a,b| bases[b][:num] <=> bases[a][:num]}
    #  alternative  = array.size - 1
    qual = bases[array[0]][:adjQual]
    oriqual = bases[array[0]][:oriQual]
    num = bases[array[0]][:num]
    refBase = @seq[ref][pos-1,1]
    @refEnv = @seq[ref][pos-7,13]
    #refEnv20=$seq[ref][pos-21,41]
    homonum = homocount(@refEnv)
    
    snpWriter.puts "#{ref}\t#{pos}\t#{refBase}\t#{homonum}\t#{@refEnv}\t#{@coverage[ref][pos]}\t#{array[0]}\t#{qual}\t#{oriqual}\t#{num}\t#{t}\t#{@snp[ref][pos]}"
  end
end
snpWriter.close

end

#Caculate the depth of coverage
def coverageDepth()
    covout = File.new(@outputFile+".ref.depth_cov", 'w')
  @coverage.keys.sort.each do |ref|
    covout.syswrite(">#{ref} depth-coverage\n")
    s = @coverage[ref].size 
    0.upto(s/80) do |i|
      covout.syswrite(@coverage[ref][(i*80+1)..(i+1)*80].join(' ') + "\n")
    end
  end
  covout.close
end
end

optHash=AtlasSubs.processArguments()
atlasSubsParser=AtlasSubs.new(optHash)
atlasSubsParser.hashNQS()
atlasSubsParser.initializeReference()
atlasSubsParser.initializeCoverage()
atlasSubsParser.parseSubs()
atlasSubsParser.printSubs()
if @coverage==true
    atlasSubsParser.coverageDepth()
end
exit(0);

