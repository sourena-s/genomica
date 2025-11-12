##This script calculates my daily commuting expenses to Donders Institute in 2025
##The input csv file (data.csv) is exported from MijnNS.NL trip history
awk -F, '
BEGIN {

IGNORECASE = 1
FPAT = "([^,]+)|(\"[^\"]+\")"
OFS = ","
  nijmegen = "(Nijmegen|Nijmegen, (Han|Radboudumc|Gymnasion|Erasmusgebouw|Spinozageb/tandheelk|Centraal Station))"
}

function connected(a, b) {
  return (a ~ nijmegen && b ~ nijmegen)
}

# Keep at most 3 Bus trips in Nijmegen per day
# one round-trip + possible one evening trip (for next day)
$10 ~ /Bus/ && connected($3, $5) {
  date = $1
  if (count[date] < 3) {
          cost = $6
  gsub(/[^0-9,]/,"",cost)
  gsub(",",".",cost)
  $6 = "€" cost
    print
    count[date]++
  }
}
' data.csv > trips.txt

awk -F, '
BEGIN {
  IGNORECASE = 1
  FPAT = "([^,]+)|(\"[^\"]+\")"
  OFS = ","
}
function to_minutes(t) {
  split(t, a, ":"); return a[1]*60 + a[2]
}
{
  date = $1
  time = $2
  from = $3
  to   = $5

# Case 1: Nijmegen -> Arnhem after 17:00 -> store temporarily
if (from ~ /Nijmegen/ && to ~ /Arnhem Centraal/ && to_minutes(time) >= 17*60) {
    evening_trip[date, time] = $0
    next
}

# Case 2: Arnhem -> Nijmegen after 17:00 -> store
if (from ~ /Arnhem Centraal/ && to ~ /Nijmegen/ && to_minutes(time) >= 17*60) {
    arnhem_trip[date, time] = $0
    next
}

# Case 3: regular trips (before 17:00) - print normally
if (to_minutes(time) < 17*60 && count[date] < 3) {
    cost = $6
    gsub(/[^0-9,]/,"",cost)
    gsub(",",".",cost)
    $6 = "€" cost
    print
    count[date]++
}

END 
{
    # Print evening Nijmegen->Arnhem trips only if no later trip exists
    for (key in evening_trip) {
        split(key, arr, SUBSEP)
        d = arr[1]; t = arr[2]

        later_trip_exists = 0
        for (akey in evening_trip) {
            split(akey, aarr, SUBSEP)
            if (aarr[1] == d && aarr[2] > t) {
                later_trip_exists = 1
                break
            }
        }

        if (!later_trip_exists) {
            print evening_trip[key]
            count[d]++
        }
    }

    # Optionally print Arnhem->Nijmegen evening trips if needed
    # for (key in arnhem_trip) print arnhem_trip[key]
}
' data.csv >> trips.txt


awk -F',' '{split($1,d,"-"); cmd="date -d \"" d[2] "/" d[1] "/" d[3] " " $2 "\" +%s"; cmd|getline e; close(cmd); print e","$0}' trips.txt |sort -k1,1n -t, > trips.sorted

awk -F',' ' BEGIN {
  IGNORECASE = 1
  FPAT = "([^,]+)|(\"[^\"]+\")"
  OFS = ","
}
{
  gsub(/€/,"",$7); gsub(/"/,"",$7); gsub(/,/, ".", $7)
  if ($2 != prev && NR>1) {
    printf "%d,Total for %s: €%.2f\n", prev_timestamp, prev, sum
    sum = 0
  }
  sum += $7
  print
  prev = $2
  prev_timestamp = $1
} END {
  if (NR>0) printf "%d,Total for %s: €%.2f\n", prev_timestamp, prev, sum
  }' trips.sorted | sort -t, -k1,1n |cut -d, -f2- > commute_data.txt

rm trips.txt trips.sorted

