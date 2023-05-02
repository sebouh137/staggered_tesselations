module plaque(lines){
    rotate(180) {
        linear_extrude(2)
            for (i=[0 : len(lines)])
            translate([10, -(i-len(lines)/2+1)*12])
                text(lines[i],size=10);
        linear_extrude(1)
            translate([100,0]) square([200,100],center=true);
    }
    
}

//plaque(["layer 0", "hole radius=XXX cm", "hole x=XXX cm"]);
//include <layer_0.scad>