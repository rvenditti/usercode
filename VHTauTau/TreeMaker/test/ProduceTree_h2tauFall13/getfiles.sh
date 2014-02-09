#!/bin/tcsh

set dirin = $*
foreach i (`ls -1 "$dirin" | grep patTuple | awk '{print $NF}'`)

echo  ''$dirin$i','

end

