foreach year ( $argv[1] $argv[2] )
        echo $year
        set j = 3
        foreach i ( 1 2 )
                echo $i
        end
        set h = `echo $i | awk '{print $1+3}'`
        echo $i,$h
end

