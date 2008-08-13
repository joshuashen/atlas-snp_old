require 'getoptlong'
require 'rsruby'

r=RSRuby.instance

opts = GetoptLong.new(
    ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
    ["--pvalue","-p", GetoptLong::REQUIRED_ARGUMENT],
    ["--help","-h", GetoptLong::REQUIRED_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

#In this freeze, the genotype determination using binomial p value is provided. For the next freeze, other ways to determine genotype will be provided. Thus, users can have multiple choices.
#Also, in next freeze, a clear log file storing the ratio of cumulative p values will be provided, users can set their own cutoff according to this information. 

if optHash.key?("--help") or !optHash.key?("--input")
  $stderr.puts "Usage: ruby __.rb -i SNP.list [ -p binomial p value cutoff ].  This program should be used after SNP evaluating."
  $stderr.puts "       -i   SNP list"
  $stderr.puts "       -p   default value: 0.001.    This option permits user to denote the cutoff for binomial p value."
  exit
end

if optHash.key?("--pvalue")
    $p=optHash["--pvalue"].to_f
else
    $p=0.001
end

file_name=File.basename(optHash["--input"])
file_path=File.dirname(optHash["--input"])

file=file_path+'/'+file_name+'.pvalue.gt'
fo=File.new(file,'w')

File.new(optHash["--input"],'r').each do |line|
    cols=line.split(/\s+/)
    
    chr=cols[0]
    refbase=cols[2].upcase
    cov=cols[5].to_i
    snpbase=cols[6].upcase
    varinum=cols[9].to_i
    othernum=cols[10].to_i
    refnum=cov-othernum
    
    coverage=refnum+varinum
    
    if chr=="refName"
        fo.print cols[0..-1].join("\t")+"\tBinomial P\tGenotype\n"
    else
    
    ratio=(refnum.to_f/coverage.to_f).to_f
    #r.eval_R("refno=refnum")
    
    h=r.eval_R("binom.test(#{refnum},#{coverage},p=0.5)")
    pvalue=h["p.value"]
    prefix=cols[0..-1].join("\t") + "\t"
    
    if refnum > 0
        if pvalue < $p
            if ratio >= 0.5
                gtrr="#{refbase}#{refbase}"
                fo.print prefix+"#{pvalue}\t#{gtrr}\n"
            else
                gt="#{snpbase}#{snpbase}"
                fo.print prefix+"#{pvalue}\t#{gt}\n"
            end
        else
            gt=[refbase, snpbase].sort.join()
            fo.print prefix+"#{pvalue}\t#{gt}\n"
        end
    else
        gt="#{snpbase}#{snpbase}"
        fo.print prefix+"#{pvalue}\t#{gt}\n"
    end
    
    end
end
fo.close


