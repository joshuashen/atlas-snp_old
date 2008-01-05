# convert the coordinates to global

$pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/

offset = 0
tname = ''

while line=ARGF.gets do
	if line.match($pattern)
		score,sub,del,ins,query,qstart,qend,qright,tstrand,target,d1,d2,d3,lab =$1.to_f,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11,$12,$13,$14
		if target =~ /^(\S+)\_(\d+)\_(\d+)/
			offset  = $2.to_i
			tname = $1
		end
		
		tstart = 0
		tend = 0
		tright = 0
      
		if d1 =~ /\((\d+)\)/
			tright, tstart, tend = $1.to_i, d2.to_i, d3.to_i
		elsif d3=~ /\((\d+)\)/
			tstart, tend, tright = d1.to_i, d2.to_i, $1.to_i
		end

		tstart = tstart + offset - 1
		tend = tend + offset - 1
	elsif	line=~ /^\s\S+/