include <segments.scad>
wall_thickness=0.8;
frame_height=3;
module connect(start, end){
    dx=start-end;
    length=norm(dx);
    phi=atan2(dx[1],dx[0]);
     translate((start+end)/2) rotate(phi) square([length, wall_thickness],true);
     $fn=12;
     translate(start) circle(r=wall_thickness/2);
     translate(start) circle(r=wall_thickness/2);
}

linear_extrude(frame_height) {
    union(){
        for(segment = segments)
            connect(segment[0],segment[1]);
        for (corner =corners)
            translate(corner) square(wall_thickness,center=true);
        }
}