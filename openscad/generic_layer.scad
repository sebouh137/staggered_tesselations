cm=10;
mm=1;

module generic_layer(z, holeX, holeR, thickness, det_height=62.24*cm, det_width=60.96*cm, det_gap=0.4*cm,$fn=48){
    translate([0,0,z]) linear_extrude(thickness,center=false) difference(){
        translate([-100, 0,0]) square([det_width, det_height], center=true);
        translate([holeX, 0,0]) circle(holeR);
        square([det_gap,det_height], center=true);
        translate([holeX/2,0,0])square([abs(holeX), holeR*2], center=true);
    }
}

//generic_layer(300, -70, 150, 16);