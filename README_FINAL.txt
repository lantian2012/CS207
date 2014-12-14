
Overview:
  Our project is completed in three phases. The result of our first phase is in the root folder. The meshed_mass_spring.cpp file demos a inflated sphere that bounces against a plane. The result of our second phase is in the folder Combine. We add the collision detector, the visualizer, and a constraint for collision between spheres in this phase.  The meshed_mass_sping.cpp demos the collision of two inflated spheres with wind. In the third phase, we further add shallow water extension. We are able to show the visual effect of a ball that drops into a pond. Files of the third phase lies in the folder Combine2.
					 				
1.tianlian-cs207: Phase 1
Feature:
  In the meshed_mass_spring.hpp file, we designed a system based on meshed mass-spring models. By using a triangular Mesh instead, we implemented forces and constraints applied to the surfaces (triangles) of the mesh. 

run the code:
  compile: make meshed_mass_sping
  run: ./meshed_mass_spring data/XXX.nodes data/XXX.tris
  For bouncing ball simulation, we suggest using sphere200.nodes and sphere200.tris for good visual effects.

interesting things you should try:
  In this folder, users can try several kinds of meshed mass-spring physics simulation:
  Flying cloth: 
    Data: users can choose data among tub2.*, tub3.*, tub4.* files.
    Forces: users can add or combine different forces: wind force, mass-spring force, gravity, and damping force.
    Constraints: users can add or combine different constraints: plane constraint, box constraint, and constant constraint.

  Bouncing ball: 
    Data: users can choose data among sphere12.*, sphere200.*, sphere1082.* files.
    Forces: users can add or combine different forces: wind force, mass-spring force, gravity, pressure force, and damping force.
    Constraints: users can add or combine different constraints: plane constraint, box constraint, and constant constraint.

2.combine: Phase 2
Feature:
  We add the collision detector, the collision constraint and a visualizer to our meshed mass spring. The demo is the collision of two balls. We detect the nodes that collide with the other mesh, find the collision plane, and then set the velocity and position of collided nodes. The two balls will bounce back after collision. As we also add wind in the demo, the two balls will tend to fly towards a direction.

run the code:
  compile: make meshed_mass_sping
  run: ./meshed_mass_spring data/sphereXXX.nodes data/sphereXXX.tris data/sphereXXX.nodes data/sphereXXX.tris
  We suggest using sphere200.nodes and sphere200.tris for good visual effects.

run-time commands or interaction:
  using keyboard to control the object movement：
    press W/S to control the movement in X-direction
    press A/D to control the movement in Y-direction
    press Z/C to control the movement in Z-direction         
  Pause/ accelerate/ decelerate/ recover the simulation by mouse and keyboard:
    press the Left Click of mouse to Pause
    press Up to accelerate
    press Down to decelerate
    press Left to recover
  Adjust the objects’ color with keyboard.
    press 1~9 to change the object’s color

interesting things you should try:
  Collide the ball for many times
    You can use the command to change the position of one ball, and make it collide with the other ball for several times. You can see that in our set-up, the air filled inside the ball is little, so colliding will cause interesting deformation of the spheres. 
  Add box constraint
    You can add the box constraint in the main function. After properly setting the position of the box, you will be able see two balls collide with each other and bounce against the wall.

3.combine2: Phase 3
		
Feature:
  In this folder, we add the shallow water extensions, the collision detector, the collision constraint and a visualizer to our meshed mass spring. The demo is the collision of an inflated ball and water surface. We detect the nodes that collide with the other mesh, find the collision plane, and then set the velocity and position of collided nodes. The inflated ball will bounce back after collision. 

run the code:
  compile: make meshed_mass_sping
  run: ./meshed_mass_spring data/sphereXXX.nodes data/sphereXXX.tris data/sphereXXX.nodes data/sphereXXX.tris
  For a good visual effects, we suggest using sphere200.nodes and sphere200.tris for the inflated ball mesh and using tub3.nodes and tub3.tris for the water surface mesh.

run-time commands or interaction:
  Pause/ accelerate/ decelerate/ recover the simulation by mouse and keyboard:
    press the Left Click of mouse to Pause
    press Up to accelerate
    press Down to decelerate
    press Left to recover
  Adjust the objects’ color with keyboard.
    press 1~9 to change the object’s color

Remark:
  Due to a defect of the Morton code module supplied by another group, this demo will reach a point when the mesh is too big and cause a failure in the Morton code module. 

