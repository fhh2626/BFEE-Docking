from __future__ import absolute_import

from openmm import *
from openmm.unit import *
# from openmm.unit import MOLAR_GAS_CONSTANT_R,angstrom

import numpy as np

class TABFIntegrator(CustomIntegrator):
    def __init__(self, dt, malpha, temperature_alpha, frictionCoeff_alpha, hill_height, hill_sigma, cutoff_width, wtm_temperature, temperature, frictionCoeff, temperature_range, bins, nfull, checkpoint_file=None):
        """Create an TABFIntegrator.

        Parameters
        ----------
        dt : time, The integration time step to use, in unit fs
        
        * Extended freedom Parameters
        malpha : mass of extended freedom, in atomic unit
        
        temperature_alpha : temperature of extended freedom, in unit Kelvin
        
        frictionCoeff_alpha : friction coefficient of extended freedom, in unit 1
        
        * Langevin Parameters
        temperature : temperature, in unit Kelvin
        
        frictionCoeff : friction coefficient, in unit 1
        
        * ABF Parameters
        temperature_range: list of temperature range
        
        bins: number of bins
        
        nfull: the minumum number of 
        
        """
        CustomIntegrator.__init__(self, dt)
        self.temperature = temperature
        self.temperature_alpha = temperature_alpha
        self.frictionCoeff_alpha = frictionCoeff_alpha
        
        # WTM parameters
        try:
            hill_height.unit.get_name()
            self.hill_height = hill_height
        except:
            self.hill_height = hill_height * kilojoule_per_mole
        self.hill_sigma = hill_sigma
        self.cutoff_width = cutoff_width
        self.wtm_temperature = wtm_temperature
        
        self.frictionCoeff = frictionCoeff
        self.temperature_range = temperature_range
        self.bins = bins
        self.nfull = nfull
        
        # processing checkpoint
        self.checkpoint_file_ = checkpoint_file
        
        if self.checkpoint_file_:
            # alpha and its velocity
            self.checkpoint_file = np.load(self.checkpoint_file_,allow_pickle=True).item()
            self.ckpt_alpha = self.checkpoint_file['ckpt_alpha']
            self.ckpt_valpha = self.checkpoint_file['ckpt_valpha']
            
            # histogram info
            self.ckpt_bin_count = self.checkpoint_file['ckpt_bin_count']
            self.ckpt_Favg = self.checkpoint_file['ckpt_Favg']
            self.ckpt_wtmpotential = self.checkpoint_file['ckpt_wtmpotential']
            self.ckpt_wtmbias = self.checkpoint_file['ckpt_wtmbias']

        # init wtmabf
        self._init_abf()
        self._init_wtm()
        
        self.addGlobalVariable("a", math.exp(-frictionCoeff*dt))
        self.addGlobalVariable("b", math.sqrt(1-math.exp(-2*frictionCoeff*dt)))
        self.addGlobalVariable("a_alpha", math.exp(-frictionCoeff_alpha*dt))
        self.addGlobalVariable("b_alpha", math.sqrt(1-math.exp(-2*frictionCoeff_alpha*dt)))
        
        from openmm.unit import MOLAR_GAS_CONSTANT_R,angstrom
        self.addGlobalVariable('kT', MOLAR_GAS_CONSTANT_R*temperature)
        self.addGlobalVariable('du', 1*angstrom)
        self.addPerDofVariable("x1", 0)
        
        self.addGlobalVariable("malpha", malpha)
        self.addGlobalVariable("Talpha", temperature_alpha)
        
        self.addGlobalVariable('k',MOLAR_GAS_CONSTANT_R)
        self.addGlobalVariable('kTalpha', MOLAR_GAS_CONSTANT_R*temperature_alpha)
        
        # initialize alpha
        if self.checkpoint_file_:
            self.addGlobalVariable("alpha", self.ckpt_alpha)
            self.addGlobalVariable("valpha",self.ckpt_valpha)

        else:
            self.addGlobalVariable("alpha", self.alpha_min)
            self.addGlobalVariable("valpha", ((MOLAR_GAS_CONSTANT_R*temperature_alpha/malpha).sqrt())*np.random.normal(size=1)/angstrom)
            
        self.addUpdateContextState()
        
        # update xyz
        self.addComputePerDof("v", "v+0.5*dt*f/(m*alpha)")
        self.addComputePerDof("x", "x+(dt/2)*v")
        self.addComputePerDof("v", "a*v + b*sqrt(kT/m)*gaussian")
        self.addComputePerDof("x", "x+(dt/2)*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions();
        self.addComputePerDof("v", "v+(x-x1)/dt");
        self.addConstrainVelocities()
        self.addComputePerDof("v", "v+0.5*dt*f/(m*alpha)")
        self.addConstrainVelocities();
        
        # update alpha
        self.__wtmabf__()
        
    def _init_wtm(self):
        # for calculating dx
        colvars = np.linspace(self.alpha_min - (1+3*self.cutoff_width)*self.alpha_window_width, self.alpha_max + 3*self.cutoff_width*self.alpha_window_width, self.bins+2+6*self.cutoff_width)
        
        dx = np.ones((self.bins+2+6*self.cutoff_width,1)) * colvars[None,...] # broadcasting
        dx = dx - colvars[..., None] # broadcasting
        dx2 = dx*dx
        self.dx = dx
        self.dx2 = dx2
        
        self.addGlobalVariable('hill_height', self.hill_height)
        self.addGlobalVariable('hill_sigma', self.hill_sigma*self.alpha_window_width)
        self.addGlobalVariable('wtm_temperature', self.wtm_temperature)
        self.addGlobalVariable('w', self.hill_height)
        
        if self.checkpoint_file_:
            for bin_id in range(self.bins+2+6*self.cutoff_width):
                self.addGlobalVariable(f"WTMpotential{bin_id-self.cutoff_width}", self.ckpt_wtmpotential[bin_id]*kilojoule_per_mole)
                self.addGlobalVariable(f"WTMbias{bin_id-self.cutoff_width}", self.ckpt_wtmbias[bin_id]*kilojoule_per_mole) # calculated by WTMpotential
            
        else:
            for bin_id in range(self.bins+2+6*self.cutoff_width):
                self.addGlobalVariable(f"WTMpotential{bin_id-self.cutoff_width}", 0*kilojoule_per_mole)
                self.addGlobalVariable(f"WTMbias{bin_id-self.cutoff_width}", 0*kilojoule_per_mole)

        self.addGlobalVariable("pot_current", 0*kilojoule_per_mole)
        self.addGlobalVariable("Fbias_wtm", 0*kilojoule_per_mole)
        self.addUpdateContextState()

    def _init_abf(self):
        alpha_min = min(self.temperature_range)/(self.temperature)
        alpha_max = max(self.temperature_range)/(self.temperature)
        alpha_window_width = (alpha_max - alpha_min)/self.bins
        
        self.alpha_min = alpha_min
        self.alpha_max = alpha_max
        self.alpha_window_width = alpha_window_width
        
        self.addGlobalVariable('alpha_min', alpha_min)
        self.addGlobalVariable('alpha_max', alpha_max)
        self.addGlobalVariable('nfull', self.nfull)
        self.addGlobalVariable('alpha_window_width', alpha_window_width)
        
        # Removing boundary effect
        self.addGlobalVariable('alpha_min_boundary', alpha_min - (1+2*self.cutoff_width)*alpha_window_width)
        self.addGlobalVariable('alpha_max_boundary', alpha_max + (1+2*self.cutoff_width)*alpha_window_width)
        
        if self.checkpoint_file_:
            for bin_id in range(self.bins+2+4*self.cutoff_width):
                self.addGlobalVariable(f"bin{bin_id}", self.ckpt_bin_count[bin_id])
                self.addGlobalVariable(f"Favg{bin_id}", self.ckpt_Favg[bin_id]*kilojoule_per_mole)
                
        else:
            for bin_id in range(self.bins+2+4*self.cutoff_width):
                self.addGlobalVariable(f"bin{bin_id}", 0)
                self.addGlobalVariable(f"Favg{bin_id}", 0*kilojoule_per_mole)
            
        self.addGlobalVariable("at_bin", 0)
        self.addGlobalVariable("Fbias", 0*kilojoule_per_mole)
        self.addGlobalVariable("Fbias_temp", 0*kilojoule_per_mole)
        self.addGlobalVariable("f_alpha", 0*kilojoule_per_mole)
        self.addUpdateContextState()
        
    def __wtmabf__(self):
        # cal bin index
        self.addComputeGlobal("f_alpha","energy/(alpha*alpha)")
        self.addComputeGlobal('at_bin','ceil((alpha - alpha_min_boundary)/alpha_window_width)')
        
        for bin_id in range(self.bins+2+4*self.cutoff_width):
            self.beginIfBlock(f"at_bin={bin_id}")
            # wtm part
            self.addComputeGlobal("w",f"hill_height*exp(-WTMpotential{bin_id}/(k*wtm_temperature))")
            
            for tmp_bin_id in range(self.cutoff_width*2+1):
                
                self.addComputeGlobal(f"pot_current",f"w*exp({-self.dx2[bin_id+self.cutoff_width,bin_id+tmp_bin_id]}/(2*hill_sigma*hill_sigma))")
                
                self.addComputeGlobal(f"WTMpotential{bin_id-self.cutoff_width+tmp_bin_id}",f"WTMpotential{bin_id-self.cutoff_width+tmp_bin_id} + pot_current")
                
                self.addComputeGlobal(f"WTMbias{bin_id-self.cutoff_width+tmp_bin_id}",f"WTMbias{bin_id-self.cutoff_width+tmp_bin_id} + pot_current*{self.dx[bin_id+self.cutoff_width,bin_id+tmp_bin_id]}/(hill_sigma*hill_sigma)")
                
            self.addComputeGlobal("Fbias_wtm",f"WTMbias{bin_id}")
            
            # abf part
            self.addComputeGlobal(f"bin{bin_id}", f"bin{bin_id} + 1")
            self.addComputeGlobal(f"Favg{bin_id}",f"(Favg{bin_id}/bin{bin_id}) * (bin{bin_id}-1) + f_alpha/bin{bin_id}")
            self.beginIfBlock(f"bin{bin_id} > nfull")
            self.addComputeGlobal("Fbias",f"-Favg{bin_id}")
            self.addComputeGlobal("Fbias_temp",f"Fbias")
            self.endBlock()
            
            self.endBlock()
            
        # updating velocity of alpha
        self.addComputeGlobal("valpha", "valpha+0.5*dt*(f_alpha+Fbias+Fbias_wtm)/(malpha*du*du)")
        self.addComputeGlobal("alpha", "alpha+(dt/2)*valpha")
        self.addComputeGlobal("valpha", "a_alpha*valpha + b_alpha*sqrt(k*Talpha/(malpha))*gaussian/du")
        self.addComputeGlobal("alpha", "alpha+(dt/2)*valpha")
        self.addComputeGlobal("f_alpha","energy/(alpha*alpha)")
        self.addComputeGlobal("valpha", "valpha+0.5*dt*(f_alpha+Fbias+Fbias_wtm)/(malpha*du*du)")
        
        self.beginIfBlock("alpha > alpha_max_boundary")
        self.addComputeGlobal("alpha","alpha_max_boundary")
#        self.addComputeGlobal("valpha", '-abs(sqrt(k*Talpha/malpha)*gaussian/du)')
        self.addComputeGlobal("valpha", '-valpha')
        self.endBlock()
        
        self.beginIfBlock("alpha < alpha_min_boundary")
        self.addComputeGlobal("alpha","alpha_min_boundary")
        self.addComputeGlobal("valpha", '-valpha')
#        self.addComputeGlobal("valpha", 'abs(sqrt(k*Talpha/malpha)*gaussian/du)')
        self.endBlock()
        
        self.addComputeGlobal("Fbias","0")
        
class CustomVariableReporter:
    def __init__(self, file, reportInterval, variable_names):
        """
        Parameters
        ----------
        file : str
            输出文件路径
        reportInterval : int
            报告间隔步数
        variable_names : str or list
            要报告的变量名称，可以是单个字符串或字符串列表
        """
        self._file = open(file, 'w')
        self._reportInterval = reportInterval
        # 确保 variable_names 是列表形式
        self._variable_names = [variable_names] if isinstance(variable_names, str) else variable_names
        
        # 写入表头
        header = "Step " + " ".join(self._variable_names)
        print(header, file=self._file)
        self._file.flush()

    def __del__(self):
        if hasattr(self, '_file'):
            self._file.close()

    def report(self, simulation, state):
        # 获取当前步数
        step = simulation.currentStep
        
        # 获取所有变量的值
        values = []
        for var_name in self._variable_names:
            value = simulation.integrator.getGlobalVariableByName(var_name)
            values.append(f"{value:.6f}")
        
        # 输出格式化结果
        output_line = f"{step} " + " ".join(values)
        print(output_line, file=self._file)
        self._file.flush()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False)
        
class PMFReporter:
    def __init__(self, file, reportInterval, bins, cutoff_width):
        self._file = file
        self._reportInterval = reportInterval
        self.bins = bins
        self.cutoff_width = cutoff_width
        self.pmf = np.tril(np.ones(shape=(bins,bins)),k=0)

    def __del__(self):  # 确保文件在结束时关闭
        pass
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        f = open(self._file, 'w')
        integrator = simulation.integrator
        f.write('alpha\tFavg\tpmf\thist\n')
        try:
            self.alpha_min
        except:
            self.alpha_min = integrator.getGlobalVariableByName('alpha_min')
            self.alpha_max = integrator.getGlobalVariableByName('alpha_max')
            self.alpha_window_width = integrator.getGlobalVariableByName('alpha_window_width')

        # cal pmf and hist 
        Favg = []
        hist = []
        for i in range(1+2*self.cutoff_width,self.bins+1+2*self.cutoff_width):
            Favg.append(integrator.getGlobalVariableByName(f'Favg{i}'))
            hist.append(integrator.getGlobalVariableByName(f'bin{i}'))
        
        # cal pmf
        F_avg = -np.array(Favg)
        pmf = self.pmf*self.alpha_window_width
        pmf = (pmf*F_avg.reshape(1,-1)).sum(1)
        pmf -= pmf.min()
        pmf *= 0.239

        for i in range(self.bins):
            print(f"%.4f\t%.4f\t%.4f\t%d"%((i+0.5)*self.alpha_window_width+self.alpha_min,-F_avg[i],pmf[i],hist[i]), file=f)

        f.flush()
        f.close()
        
class MTDReporter:
    def __init__(self, file, reportInterval, bins, cutoff_width):
        self._file = file
        self._reportInterval = reportInterval
        self.bins = bins
        self.cutoff_width = cutoff_width

    def __del__(self):  # 确保文件在结束时关闭
        pass
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        f = open(self._file, 'w')
        integrator = simulation.integrator
        f.write('alpha\tWTMpotential\tWTMbias\n')
        try:
            self.alpha_min
        except:
            self.alpha_min = integrator.getGlobalVariableByName('alpha_min')
            self.alpha_max = integrator.getGlobalVariableByName('alpha_max')
            self.alpha_window_width = integrator.getGlobalVariableByName('alpha_window_width')

        WTMpotential = []
        WTMbias = []
        for bin_id in range(2*self.cutoff_width+1,self.bins+2*self.cutoff_width+1):
            WTMpotential.append(integrator.getGlobalVariableByName(f"WTMpotential{bin_id}"))
            WTMbias.append(integrator.getGlobalVariableByName(f"WTMbias{bin_id}"))

        for i in range(self.bins):
            print(f"%.4f\t%.4f\t%.4f"%((i+0.5)*self.alpha_window_width+self.alpha_min,WTMpotential[i],WTMbias[i]), file=f)

        f.flush()
        f.close()

class CKPTReporter:
    def __init__(self, file, reportInterval, bins, cutoff_width):
        self._file = file
        self._reportInterval = reportInterval
        self.bins = bins
        self.cutoff_width = cutoff_width

    def __del__(self):  # 确保文件在结束时关闭
        pass
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        ckpt_dictionary = {}
        integrator = simulation.integrator
        
        ckpt_dictionary['ckpt_step'] = simulation.currentStep
        ckpt_dictionary['ckpt_alpha'] = integrator.getGlobalVariableByName('alpha')
        ckpt_dictionary['ckpt_valpha'] = integrator.getGlobalVariableByName('valpha')
        ckpt_dictionary['ckpt_bin_count'] = []
        ckpt_dictionary['ckpt_Favg'] = []
        ckpt_dictionary['ckpt_wtmpotential'] = []
        ckpt_dictionary['ckpt_wtmbias'] = []
        
        for bin_id in range(self.bins+2+4*self.cutoff_width):
            ckpt_dictionary['ckpt_bin_count'].append(integrator.getGlobalVariableByName(f"bin{bin_id}"))
            ckpt_dictionary['ckpt_Favg'].append(integrator.getGlobalVariableByName(f"Favg{bin_id}"))
            
        for bin_id in range(self.bins+2+6*self.cutoff_width):   
            ckpt_dictionary['ckpt_wtmpotential'].append(integrator.getGlobalVariableByName(f"WTMpotential{bin_id-self.cutoff_width}"))
            ckpt_dictionary['ckpt_wtmbias'].append(integrator.getGlobalVariableByName(f"WTMbias{bin_id-self.cutoff_width}"))
       
        np.save(self._file, ckpt_dictionary)