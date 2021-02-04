import bpy
import glob

if __name__ == "__main__":
    for i in range(0,1000,4):
        bpy.ops.object.select_all(action='SELECT')
        bpy.ops.object.delete()
        files = glob.glob('/home/backes/dev/pbs20-wall/recordings/rec3/obj_*_'+str(i)+'.obj')
        for file in files:
            if (str(file) == '/home/backes/dev/pbs20-wall/recordings/rec3/obj_74_step_'+str(i)+'.obj') :
                continue
            if (str(file) == '/home/backes/dev/pbs20-wall/recordings/rec3/obj_75_step_'+str(i)+'.obj') :
                continue
            bpy.ops.import_scene.obj(filepath=str(file))
        bpy.ops.export_scene.obj(filepath='/home/backes/dev/pbs20-wall/recordings/out_1/wall_frame_'+str(i)+'.obj', use_materials=False)
        bpy.ops.object.select_all(action='SELECT')
        bpy.ops.object.delete()
