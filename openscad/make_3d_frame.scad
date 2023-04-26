module connect(start, end, thickness=0.41){
    dx=start-end;
    length=norm(dx);
    phi=atan2(dx[1],dx[0]);
     translate((start+end)/2) rotate(phi) square([length, thickness],true);
     $fn=12;
     translate(start) circle(r=thickness/2);
     translate(start) circle(r=thickness/2);
    
}

//segments = [
//    [[0,0],[0,1]],
//    [[0,0.5],[1,0.5]]
//    ];
//segments_file="segments.scad";
include <tmp.scad>
linear_extrude(3) {
    union(){
        for(segment = segments)
            connect(segment[0],segment[1]);
        }
}
