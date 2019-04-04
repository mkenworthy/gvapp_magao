def cen_rot (im, rot, dim_out, offset1=(0,0), offset2=(0,0), order=2):
    """cen_rot - takes a cube of images im, and a set of rotation angles in rot,
                and translates the middle of the frame with a size dim_out to the middle of
                a new output frame with an additional rotation of rot.
                """
    import numpy as np
    from scipy.ndimage import affine_transform
    a = rot * np.pi / 180.
    
    # make a rotation matrix
    transform=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
    
    # calculate total offset for image output

    # determine centre of input image

    # -0.5 is there for fencepost counting error
    c_in = np.array(offset1) - 0.5
    c_out = 0.5 * np.array(dim_out) - 0.5

    # c_out has to be pre-rotated to make offset correct
    offset = c_in - c_out.dot(transform) - np.array(offset2).dot(transform)
    
    # perform the transformation
    dst=affine_transform( \
        im,transform.T, order=order,offset=offset, \
        output_shape=dim_out, cval=0.0)
    return(dst)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.visualization import quantity_support

    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    t1 = np.random.random((11,17))
    #t1 = np.ones((11,17))
    ax1.imshow(t1)

    nx, ny = t1.shape
    t2 = cen_rot(t1,30,(25,29), offset1=(0.5*nx, 0.5 * ny))
    ax2.imshow(t2)
    plt.show()
    plt.draw()
