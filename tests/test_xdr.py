#!/usr/bin/env python
import pytest
import numpy as np
from pmx.model import Model
from pmx.xtc import Trajectory
from numpy.testing import assert_almost_equal

@pytest.mark.parametrize("struct, traj", [
    ('peptide.pdb', 'peptide.trr'),
    ('peptide.pdb', 'peptide.xtc')
])
def test_trajectory_xdr(gf, struct, traj):
    m = Model(gf(struct))
    t = Trajectory(gf(traj))
    for frame in t:
        frame.update(m)
    assert_almost_equal(m.atoms[2].x[1], 34.70, decimal=2)
    t.close()
    
@pytest.mark.parametrize("struct, traj", [
    ('peptide_2.pdb', 'peptide_2.trr'),
])    
def test_trajectory_xdr_w_vel_read_write(gf, struct, traj, tmpdir):
    to_file=str(tmpdir)+"/test.trr"
    m = Model(gf(struct))
    t = Trajectory(gf(traj))
    to = Trajectory(to_file, mode='Out', atomNum = len(m.atoms))
    
    #copy trajectory
    for frame in t:
        frame.update(m, uv=True, uf=True)
        
        x = np.zeros(len(m.atoms)*3)
        v = np.zeros(len(m.atoms)*3)
        f = np.zeros(len(m.atoms)*3)
        for i, atom in enumerate(m.atoms):
            x[i*3:(i+1)*3]=atom.x
            v[i*3:(i+1)*3]=atom.v
            f[i*3:(i+1)*3]=atom.f
        to.write_xtc_frame(step=frame.step, time=frame.time,
                                        lam=1.0, box=frame.box, x=x, v=v, f=f,
                                        units=m.unity, bTrr=True )
        
    #    1GLU     H2    3  10.955   8.569   4.715 -2.4632  0.5158  1.8410
    assert_almost_equal(m.atoms[2].x[1], 85.69, decimal=2)
    assert_almost_equal(m.atoms[2].v[1], 5.158, decimal=2)
    t.close()
    to.close()
    
    #read written trajectory
    t = Trajectory(to_file)
    for frame in t:
        frame.update(m, uv=True, uf=True)
        
    #    1GLU     H2    3  10.955   8.569   4.715 -2.4632  0.5158  1.8410
    assert_almost_equal(m.atoms[2].x[0], 109.55, decimal=2)
    assert_almost_equal(m.atoms[2].v[0], -24.632, decimal=2)
    t.close()