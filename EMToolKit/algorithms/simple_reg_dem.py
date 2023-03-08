# +
# Just repeating the IDL comment block here for now:
# [dems,chi2] = simple_reg_dem(data, errors, exptimes, logt, tresps,
# kmax=100, kcon=5, steps=[0.1,0.5], drv_con=8.0, chi2_th=1.0, tol=0.1):
#
# Computes a DEM given a set of input data, errors, exposure times, temperatures, and
# temperature response functions. Inputs:
# data: image data (dimensions n_x by n_y by n_channels)
# errors: image of uncertainties in data (dimensions n_x by n_y by n_channels)
# exptimes: exposure times for each image (dimensions n_channels)
# logt: Array of temperatures (dimensions n_temperatures)
# tresps: Arrays of temperature response functions (dimensions n_temperatures by n_channels)
#
# Returns:
# array of DEMs (dimensions n_x by n_y by n_temperatures). Dimensions depend on the units
# of logt (need not be log temperature, despite the name) and tresps. For instance,
# if tresps has units cm^5 and the logt values are base 10 log of temperature in Kelvin,
# the output DEM will be in units of cm^{-5} per unit Log_{10}(T).
#
# Outputs:
# chi2: Array of reduced chi squared values (dimensions n_x by n_y)
#
# Optional Keywords:
# kmax: The maximum number of iteration steps (default - 100).
# kcon: Initial number of steps before terminating if chi^2 never improves (default - 5).
# steps: Two element array containing large and small step sizes (default - 0.1 and 0.5).
# drv_con: The size of the derivative constraint - threshold delta_0 limiting the change in
# log DEM per unit input temperature (default - 8; e.g., per unit log_{10}(T)).
# chi2_th: Reduced chi^2 threshold for termination (default - 1).
# tol: How close to chi2_th the reduced chi^2 needs to be (default - 0.1).
#
# Author: Joseph Plowman -- 09-15-2021
# See: https://ui.adsabs.harvard.edu/abs/2020ApJ...905...17P/abstract
# -
def simple_reg_dem(data, errors, exptimes, logt, tresps,
                   kmax=100, kcon=5, steps=[0.1, 0.5], drv_con=8.0, chi2_th=1.0, tol=0.1):

    import numpy as np
    from scipy.linalg import cho_factor, cho_solve

    [nt, nd] = tresps.shape
    nt_ones = np.ones(nt)
    [nx, ny, nd] = data.shape
    dT = logt[1:nt]-logt[0:nt-1]
    [dTleft, dTright] = [
        np.diag(np.hstack([dT, 0])), np.diag(np.hstack([0, dT]))]
    [idTleft, idTright] = [
        np.diag(np.hstack([1.0/dT, 0])), np.diag(np.hstack([0, 1.0/dT]))]
    Bij = ((dTleft+dTright)*2.0 + np.roll(dTright, -
           1, axis=0) + np.roll(dTleft, 1, axis=0))/6.0
    Rij = np.matmul((tresps*np.outer(nt_ones, exptimes)).T,
                    Bij)  # Matrix mapping coefficents to data
    Dij = idTleft+idTright - \
        np.roll(idTright, -1, axis=0) - np.roll(idTleft, 1, axis=0)
    regmat = Dij*nd/(drv_con**2*(logt[nt-1]-logt[0]))
    rvec = np.sum(Rij, axis=1)

    dems = np.zeros([nx, ny, nt])
    chi2 = np.zeros([nx, ny]) - 1.0
    for i in range(0, nx):
        for j in range(0, ny):
            err = errors[i, j, :]
            dat0 = np.clip(data[i, j, :], 0.0, None)
            s = np.log(np.sum((rvec)*(np.where(dat0 < 1.0e-2, 1.0e-2, dat0)/err**2)) /
                       np.sum((rvec/err)**2)/nt_ones)
            for k in range(0, kmax):
                # Correct data by f(s)-s*f'(s)...
                dat = (dat0-np.matmul(Rij, ((1-s)*np.exp(s))))/err
                mmat = Rij*np.outer(1.0/err, np.exp(s))
                amat = np.matmul(mmat.T, mmat)+regmat
                try:
                    [c, low] = cho_factor(amat)
                except:
                    break
                c2p = np.mean((dat0-np.dot(Rij, np.exp(s)))**2/err**2)
                deltas = cho_solve((c, low), np.dot(mmat.T, dat))-s
                deltas *= np.clip(np.max(np.abs(deltas)), None,
                                  0.5/steps[0])/np.max(np.abs(deltas))
                # Direction sign; is chi squared too large or too small?
                ds = 1-2*(c2p < chi2_th)
                c20 = np.mean(
                    (dat0-np.dot(Rij, np.exp(s+deltas*ds*steps[0])))**2.0/err**2.0)
                c21 = np.mean(
                    (dat0-np.dot(Rij, np.exp(s+deltas*ds*steps[1])))**2.0/err**2.0)
                interp_step = (
                    (steps[0]*(c21-chi2_th)+steps[1]*(chi2_th-c20))/(c21-c20))
                s += deltas*ds*np.clip(interp_step, steps[0], steps[1])
                chi2[i, j] = np.mean((dat0-np.dot(Rij, np.exp(s)))**2/err**2)
                if ((ds*(c2p-c20)/steps[0] < tol)*(k > kcon) or np.abs(chi2[i, j]-chi2_th) < tol):
                    break
            dems[i, j, :] = np.exp(s)

    return dems, chi2


def simple_reg_dem_cuda(data, errors, exptimes, logt, tresps,
                        kmax=100, kcon=5, steps=[0.1, 0.5], drv_con=8.0, chi2_th=1.0, tol=0.1):

    import numpy as np
    import torch
    import torch.nn.functional as F

    dtype = torch.float64

    data = torch.as_tensor(data).cuda().to(dtype)
    device = data.device
    errors = torch.as_tensor(errors).cuda().to(dtype)
    logt = torch.as_tensor(logt).cuda().to(dtype)
    exptimes = torch.as_tensor(exptimes).cuda().to(dtype)
    tresps = torch.as_tensor(tresps).cuda().to(dtype)

    nt = len(logt)
    [nx, ny, nd] = data.shape

    dT = logt[1:] - logt[:-1]

    # units: logK
    # triangular basis functions for integration
    Bij = (torch.diag(F.pad(dT, (1, 0)) + F.pad(dT, (0, 1))) * 2.0
           + torch.roll(torch.diag(F.pad(dT, (1, 0))), -1)
           + torch.roll(torch.diag(F.pad(dT, (0, 1))), 1)) / 6.0

    # units: DN cm^5 pix^-1 logK
    # this matrix is used to integrate the DEM to produce the model image
    # Matrix mapping coefficients to data
    Rij = (tresps * (torch.ones(nt, device=device, dtype=dtype)[:, None] @ exptimes[None, :])).T @ Bij

    Dij = torch.diag(F.pad(1 / dT, (1, 0)) + F.pad(1 / dT, (0, 1))) \
          - torch.roll(torch.diag(F.pad(1 / dT, (1, 0))), -1) \
          - torch.roll(torch.diag(F.pad(1 / dT, (0, 1))), 1)

    regmat = Dij * len(exptimes) / (drv_con ** 2 * (logt[-1] - logt[0]))

    rvec = torch.sum(Rij, axis=-1)[None, :]


    errors = errors.reshape(ny * nx, nd)
    data = data.reshape(ny * nx, nd)
    output = torch.zeros((nx * ny, nt), device=device, dtype=dtype)
    chi2 = torch.zeros((nx * ny), device=device, dtype=dtype)

    batch_size = int(2 ** 15)  # about 32k

    for i_batch in range(int(np.ceil(len(data) / batch_size))):
        errors_batch = errors[i_batch * batch_size:(i_batch + 1) * batch_size]
        data_batch = data[i_batch * batch_size:(i_batch + 1) * batch_size]
        data_batch = torch.where(data_batch < 0, 0, data_batch)
        mask = torch.ones(data_batch.shape[0], device=device, dtype=bool)

        s = torch.log(
                torch.sum(rvec * (torch.where(data_batch < 1.0e-2, 1.0e-2, data_batch) / errors_batch ** 2), dim=-1) \
                    / torch.sum((rvec / errors_batch) ** 2, dim=-1)
            ).unsqueeze(1).expand(-1, nt)

        output_batch = torch.exp(s)
        updated_model = (Rij[None, :, :] @ output_batch[:, :, None]).squeeze(-1)

        chi2_batch = torch.zeros(output_batch.shape[0], device=device, dtype=dtype) - 1.0

        indices = torch.arange(len(output_batch), device=device)

        for k in range(kmax):
            # Correct data_batch by f(s)-s *f'(s)...
            dat = (data_batch - (Rij @ ((1 - s) * torch.exp(s))[:, :, None]).squeeze(-1)) / errors_batch
            # Weight mapping by 1/err and f'(s)...
            mmat = Rij * ((1 / errors_batch)[:, :, None] @ torch.exp(s)[:, None, :])
            cholesky = torch.linalg.cholesky(mmat.transpose(1, 2) @ mmat + regmat)
            solution = torch.cholesky_solve((mmat.transpose(1, 2) @ dat[:, :, None]), cholesky).squeeze(-1)
            c2p = torch.mean((data_batch - updated_model) ** 2.0 / errors_batch ** 2, axis=-1)

            deltas = solution - s
            deltas *= (torch.clip(torch.amax(torch.abs(deltas), axis=-1), max=0.5 / steps[0]) / torch.amax(torch.abs(deltas), axis=-1))[:, None]
            ds = torch.where(c2p < chi2_th, -1, 1)  # Direction sign; is chi∧2 too large or too small?

            model0 = (Rij @ (torch.exp(s + deltas * ds[:, None] * steps[0]))[:, :, None]).squeeze(-1)
            model2 = (Rij @ (torch.exp(s + deltas * ds[:, None] * steps[1]))[:, :, None]).squeeze(-1)

            c20 = torch.mean((data_batch - model0) ** 2 / errors_batch ** 2, axis=-1)
            c21 = torch.mean((data_batch - model2) ** 2 / errors_batch ** 2, axis=-1)

            interp_step = ((steps[0] * (c21 - chi2_th) + steps[1] * (chi2_th - c20)) / (c21 - c20))

            s = s + deltas * (ds * torch.clip(interp_step, steps[0], steps[1]))[:, None]

            dems = torch.exp(s)
            updated_model = (Rij @ dems[:, :, None]).squeeze(-1)
            chi2_tmp = torch.mean((data_batch - updated_model) ** 2.0 / errors_batch ** 2, axis=-1)

            # condition 1: improvement
            # condition 2: not closer to target than tol
            new_mask = ~(((ds * (c2p - c20) / steps[0]) < tol) & (k > kcon)) \
                    & ~(torch.abs(chi2_tmp - chi2_th) < tol)

            if not torch.any(new_mask):
                break

            mask[mask.clone()] = new_mask
            data_batch = data_batch[new_mask]
            updated_model = updated_model[new_mask]
            s = s[new_mask]
            errors_batch = errors_batch[new_mask]
            indices = indices[new_mask]
            output_batch[indices] = dems[new_mask]
            chi2_batch[indices] = chi2_tmp[new_mask]

        output[i_batch * batch_size:(i_batch + 1) * batch_size] = output_batch
        chi2[i_batch * batch_size:(i_batch + 1) * batch_size] = chi2_batch
        
    output = output.reshape(nx, ny, nt).cpu().numpy()
    chi2 = chi2.reshape(nx, ny).cpu().numpy()

    return output, chi2