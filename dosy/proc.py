# SyracuseStuff/SyracuseStuff/Python Work/Franck Collab/DOSY/190104_CdSe_QD_DOSY.ipynb
import typing as T


def get_grad_list(
    n_indirect: int,
    scale_by: float = 0.535,
) -> T.List[float]:
    return (linspace(0.02, 0.98, n_indirect)) * scale_by


def get_n_indirect(data: T.Any) -> int:
    return data.get_prop('acq')['L'][22] 


def get_swh(data: T.Any)
    sw = data.get_prop('acq')['SW']
    sfo1 = data.get_prop('acq')['SFO1']
    return sw * sfo1


def load_and_init_file(
    experiment_name: str,
    experiment_type: str,
    experiment_number: int,
    N: int = 40,
    dw: int = 45,
):
    data = find_file(
        experiment_name,
        exp_type=experiment_type,
        expno=experiment_number,
    )
    n_indirect = get_n_indirect(data)

    SWH = get_swh(data)

    data.setaxis('indirect', None)
    data.chunk(
        'indirect',
        ['indirect', 'ph1&ph3_1', 'ph1&ph3_2', 'ph4'],
        [n_indirect, 2, 2, 2],
    )
    data.setaxis('ph1&ph3_1', r_[1., 2.] / 4)
    data.setaxis('ph1&ph3_2', r_[2., 3.] / 4)
    data.setaxis('ph4', r_[0, 2.] / 4)

    grad_list = get_grad_list(n_indirect)
    data.setaxis('indirect', grad_list)

    data.ift(['ph1&ph3_1', 'ph1&ph3_2', 'ph4'])
    data = data['ph1&ph3_1', 1]['ph1&ph3_2', 1]['ph4', 0].C
    data.ft('t2', shift=True)

    s = data['indirect',0].C
    x_0 = nddata(r_[-0.5:0.5:N * 1j], 'phi0')
    x_0.set_units('phi0', 'cyc')
    phi0 = exp(1j * 2 * pi * x_0)
    x_1 = (
        nddata(
            r_[-5e-2 * dw / SWH / 2:5.0e-2 * dw / SWH / 2: N * 1j],
            'phi1',
        )
        .set_units('phi1', 's')
    )
    phi1 = 1j * 2 * pi * x_1
    time = s.fromaxis('t2')
    time_new = time.data
    time_new = time_new.reshape(1, 1, time_new.shape[0])
    phi0_new = phi0.data
    phi0_new = phi0_new.reshape(phi0_new.shape[0], 1, 1)
    phi1_new = phi1.data
    phi1_new = phi1_new.reshape(1, phi1_new.shape[0], 1)
    s_new = s.data
    s_new = s_new.reshape(1, 1, s_new.shape[0])
    S = s_new * phi0_new
    s_exp = phi1_new * time_new
    S = S * (exp(s_exp))
    d_absr = S
    d_absr = abs(d_absr.real)
    d_absr = d_absr.sum(axis=2)
    ndCost = nddata(d_absr, d_absr.shape, ['phi0','phi1'])
    ndCost.setaxis('phi0', x_0.getaxis('phi0') * (1e3))
    ndCost.setaxis('phi1', x_1.getaxis('phi1') * (1e6))
