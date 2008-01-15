## input: cross_match -discrep_lists output

# note: this is for one-block alignment only
#     ignoring 1bp indels

require 'getoptlong'

opts = GetoptLong.new(
    ["--crossmatch", "-x", GetoptLong::REQUIRED_ARGUMENT],

    ["--minscore", "-m", GetoptLong::OPTIONAL_ARGUMENT],
    ["--maxsub", "-s", GetoptLong::OPTIONAL_ARGUMENT],
		["--maxindel", "-g", GetoptLong::OPTIONAL_ARGUMENT],
		["--cutoff", "-c", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end



class Indel 
  attr_accessor :ref, :query, :tstart, :tend, :qstart, :qend, :dir,:sub, :dels,:ins,:envs, :type
  
  def initialize(query, ref, qstart, qend, tstart, tend, dir,sub,dels,ins,l,envs,type)
    @query,@ref,@qstart,@qend, @tstart, @tend, @dir,@sub,@dels,@ins,@length,@envs,@type = query, ref, qstart, qend, tstart, tend, dir,sub,dels,ins,l,envs,type
  end

  def dump
    ts = @tend - @tstart 
    qs = @qend - @qstart
    gap = ts
    if @type == 'D'
      gap = ts
    else
      gap  = qs
    end
    puts "#{type}\t#{gap}\t#{@ref}\t#{@tstart}\t#{@tend}\t#{@dir}\t#{@query}\t#{@qstart}\t#{@qend}\t#{@sub}\t#{@dels}\t#{@ins}\t#{@length}\t#{envs}"
  end
end

$indels = []

pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/

sub,query, target, dir = 0.0,'','',0
dels,ins = 0.0, 0.0
length = 0
while line=ARGF.gets do
  if line.match(pattern)
    score,sub,dels,ins,query,qstart,qend,qright,tstrand,target,d1,d2,d3,lab =$1.to_f,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11,$12,$13,$14
    
    if d1 =~ /\((\d+)\)/
      tright, tstart, tend = $1.to_i, d2.to_i, d3.to_i
    elsif d3=~ /\((\d+)\)/
      tstart, tend, tright = d1.to_i, d2.to_i, $1.to_i
    end
      
    length = qright + qend
    if tstart < tend
      dir = 1
    else
      dir = -1
      tstart, tend = tend, tstart  # wicked 
    end
  elsif line=~ /^\s(\S+)\s+(\d+)\s+(\S+)\(\d+\)\s+(\d+)\s+(\S+)/
    type, qplace, tplace,envs = $1, $2.to_i, $4.to_i,$5
#    if type =~ /^D\-(\d+)/
    if type =~ /^D/
      # deletion on read
      qstart = qplace 
      qend = qstart 
      block = 1
      if type =~ /^D\-(\d+)/
        block = $1.to_i
      end
      tstart = tplace 

      if dir < 0  #
        tstart = tstart - block + 1
      end
      tend = tstart + block

      indel = Indel.new(query,target, qstart,qend, tstart, tend, dir,sub,dels,ins,length,envs, "D")
      indel.dump
     # $indels << indel
    elsif type =~ /^I/
      qstart = qplace 
      tstart = tplace  
      block = 1
      if type =~ /^I\-(\d+)/
        block = $1.to_i
      end
      qend = qstart + block
     
      if dir < 0
        tstart -= 1
      end
      tend = tstart
      indel = Indel.new(query,target, qstart,qend, tstart, tend, dir,sub,dels,ins,length,envs,"I")
      indel.dump
    end
  end
end

# $indels.each {|i| i.dump}
