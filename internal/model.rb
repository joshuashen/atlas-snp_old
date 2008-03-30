## try to model the accurary of SNPs reported by single reads

snplist = ARGV[0]
# xm = ARGV[1]

transition = {'C' => 'T', 'A' => 'G', 'T' => 'C', "G" => 'A'} 

snps = {}

num = 0


# count the occurence of homo-polymer runs in a short sequence
def homocount(str)
  s = str.size
  lasts = ''
  homo = {}
  st = 0
  flag = 0
  homostart = 0
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
      homostart = i
    end
    lasts = cur
  end
  if homo.size > 0
    top = homo.keys.sort {|a,b| homo[b] <=> homo[a]}[0]
    xx = homo[top]
  else
    xx = 1
  end
  return [xx, homostart]
end


# puts "order\tbaseChangeType\t"

puts "order\tswap\tbaseChangeType\tHomoPolymerSize\tHomoPonRead\tEnv\tQual\tDist3\tTail\tSub\tIndel\tScoreRatio\tReadSize"
File.new(snplist, 'r').each do |line|
  cols = line.split(/\t/)
  
  refbase,homo, env, snpbase, info = cols[2],cols[3].to_i, cols[4],cols[6],cols[-1]
  
#  $stderr.puts "#{refbase}\t#{snpbase}"
  
  if snpbase == transition[refbase]
    trans = 1
  else
    trans = -1  # transversion
  end

  homo,homostart = homocount(env)


  info.split(';').each do |r|
    if r =~ /^(\S)\((\d+)\)(\S+)\((\d+)\)\((\S+)\/(\d+)\)([+|-])(\S+)\((\S+)\/(\S+)\/(\d+)\)(\S+)/
      qual,name,dist,score,size, dir, renv,subs,indels,tail, kind = $2.to_i, $3, $4.to_i, $5.to_f, $6.to_i, $7, $8, $9.to_f,$10.to_f, $11.to_i, $12

      dist5 = size - dist
      rhomo, rhomos = homocount(renv)      
      if kind== 'swap'
        swap = 1
      elsif kind == 'mnp'
        swap = 0.5
      else
        swap = 0
      end
      snp = {:baseChange => trans, :homoSize => homo, :env => env, :qual => qual, :readName => name,  :dist => dist, :readScore => score,  :readSize => size, :subs => subs, :indels => indels, :tail =>tail, :sawp =>  swap }
      
      if dir == '-'
        homostart = 12 - homostart
      end

      if homostart > 6 
        homol = 0
        homor = homo
      elsif homostart + homo < 6
        homol = homo
        homor = 0
      end

      if homostart < 6 and homostart + homo > 6
        homoc = homo
      else
        homoc = 0
      end


#       puts "#{num}\t#{trans}\t#{homol}\t#{homoc}\t#{homor}\t#{env}\t#{qual}\t#{dist}\t#{score/size}\t#{size}"
      puts "#{num}\t#{swap}\t#{trans}\t#{homo}\t#{rhomo}\t#{env}\t#{qual}\t#{dist}\t#{tail}\t#{subs}\t#{indels}\t#{score/size}\t#{size}"
      
    end
  end
  num += 1
end


