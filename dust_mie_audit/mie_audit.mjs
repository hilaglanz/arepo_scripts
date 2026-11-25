// Offline audit only; no AREPO code or snapshot changes.
// Homogeneous-sphere Mie series, downward logarithmic derivative and Miller psi.
// Reference algorithms/tests: https://github.com/scottprahl/miepython
// Adapted recurrence algorithms are covered by the accompanying MIEPYTHON_LICENSE.txt.
// Negative imaginary refractive index convention; double precision throughout.
import fs from 'node:fs/promises';
const add=(a,b)=>[a[0]+b[0],a[1]+b[1]];
const sub=(a,b)=>[a[0]-b[0],a[1]-b[1]];
const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
const scale=(a,s)=>[a[0]*s,a[1]*s];
const abs2=a=>a[0]*a[0]+a[1]*a[1];
const div=(a,b)=>scale([a[0]*b[0]+a[1]*b[1],a[1]*b[0]-a[0]*b[1]],1/abs2(b));
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1];
const inv=a=>div([1,0],a);

function logarithmicDerivative(z,order) {
    let twiceInverse=div([2,0],z);
    const first=scale(twiceInverse,order+0.5);
    let term=scale(twiceInverse,-order-1.5);
    let numerator=add(term,inv(first)), denominator=term;
    let ratio=div(numerator,denominator), product=mul(first,ratio);
    let iterations=0;
    while(Math.abs(Math.sqrt(abs2(ratio))-1)>1e-12) {
        if(++iterations>10000) throw Error('Continued fraction did not converge');
        term=sub(twiceInverse,term);
        numerator=add(inv(numerator),term);
        denominator=add(inv(denominator),term);
        ratio=div(numerator,denominator);
        twiceInverse=scale(twiceInverse,-1);
        product=mul(ratio,product);
    }
    return sub(product,div([order,0],z));
}

export function mie(refractiveIndex,sizeParameter) {
    const x=sizeParameter,m=refractiveIndex;
    if(!(x>0 && m[0]>0 && m[1]<=0)) throw Error('Require x>0, n>0, k>=0');
    const terms=Math.floor(x+4.05*Math.pow(x,0.33333)+2);
    const start=terms+25+Math.floor(1.5*Math.sqrt(x));
    const psi=new Float64Array(start+2);
    psi[start]=1e-50;
    for(let n=start;n>0;n--) {
        psi[n-1]=(2*n+1)/x*psi[n]-psi[n+1];
        if(Math.abs(psi[n-1])>1e200)
            for(let k=n-1;k<psi.length;k++) psi[k]/=1e200;
    }
    const exactPsi0=Math.sin(x),exactPsi1=Math.sin(x)/x-Math.cos(x);
    const normalization=Math.abs(exactPsi1)>Math.abs(exactPsi0)
        ?exactPsi1/psi[1]:exactPsi0/psi[0];
    for(let n=0;n<=terms;n++) psi[n]*=normalization;
    const z=scale(m,x),derivatives=[];
    let derivative=logarithmicDerivative(z,terms+2);
    for(let n=terms+2;n>0;n--) {
        const nOverZ=div([n,0],z);
        derivative=sub(nOverZ,inv(add(derivative,nOverZ)));
        derivatives[n-1]=derivative;
    }
    let chiPrevious=Math.cos(x),chi=Math.cos(x)/x+Math.sin(x);
    const electric=[],magnetic=[];
    for(let n=1;n<=terms;n++) {
        const outgoingPrevious=[psi[n-1],chiPrevious];
        const outgoing=[psi[n],chi];
        const electricTerm=add(div(derivatives[n],m),[n/x,0]);
        const magneticTerm=add(mul(derivatives[n],m),[n/x,0]);
        electric.push(div(sub(scale(electricTerm,psi[n]),[psi[n-1],0]),
                          sub(mul(electricTerm,outgoing),outgoingPrevious)));
        magnetic.push(div(sub(scale(magneticTerm,psi[n]),[psi[n-1],0]),
                          sub(mul(magneticTerm,outgoing),outgoingPrevious)));
        [chiPrevious,chi]=[chi,(2*n+1)/x*chi-chiPrevious];
    }
    let extinction=0,scattering=0,angularMoment=0,back=[0,0];
    for(let j=0;j<terms;j++) {
        const n=j+1,a=electric[j],b=magnetic[j];
        extinction+=(2*n+1)*(a[0]+b[0]);
        scattering+=(2*n+1)*(abs2(a)+abs2(b));
        angularMoment+=(2*n+1)/(n*(n+1))*dot(a,b);
        if(j+1<terms) angularMoment+=n*(n+2)/(n+1)*
            (dot(a,electric[j+1])+dot(b,magnetic[j+1]));
        back=add(back,scale(sub(a,b),(2*n+1)*(n%2?-1:1)));
    }
    extinction*=2/(x*x); scattering*=2/(x*x);
    const asymmetry=scattering>0?4*angularMoment/(x*x*scattering):0;
    return {extinction,scattering,absorption:extinction-scattering,
            asymmetry,backscatter:abs2(back)/(x*x)};
}

export function validateMie() {
    const cases=[
        [[1.5,-0.5],2.5,[2.562873497454734,1.097071819088392,0.748905978948507],1e-7],
        [[1.55,-0.1],2*Math.PI*0.525/0.6328,[2.86165188243,1.66424911991,0.80128972639],1e-7],
        [[1.55,0],2*Math.PI*0.525/0.6328,[3.10543,3.10543,0.63314],1e-5],
        [[1.5,-1],1,[2.336321,0.6634538,0.192136],1e-6],
        [[1.33,-1e-5],1,[null,0.093923,0.184517],1e-6],
        [[1.33,-1e-5],100,[null,2.096594,0.868959],1e-6],
        [[1.5,-0.1],0.0665,[0.0132877,0.0000047,0.0008752],1e-7]
    ];
    const results=cases.map(([m,x,expected,tolerance])=>{
        const result=mie(m,x);
        const actual=[result.extinction,result.scattering,result.asymmetry];
        const errors=actual.map((v,i)=>expected[i]===null?0:Math.abs(v-expected[i]));
        if(Math.max(...errors)>tolerance) throw Error(JSON.stringify({m,x,result,expected,errors}));
        return {m,x,result,maximumAbsoluteError:Math.max(...errors),tolerance};
    });
    const tinyX=1e-4,m=[1.5,-0.1],mSquared=mul(m,m);
    const polarizability=div(sub(mSquared,[1,0]),add(mSquared,[2,0]));
    const small=mie(m,tinyX);
    const rayleighScattering=8/3*tinyX**4*abs2(polarizability);
    const rayleighAbsorption=-4*tinyX*polarizability[1];
    if(Math.abs(small.scattering/rayleighScattering-1)>1e-6 ||
       Math.abs(small.absorption/rayleighAbsorption-1)>1e-6) throw Error('Rayleigh check failed');
    return {referenceCases:results,rayleigh:{small,rayleighScattering,rayleighAbsorption}};
}

export async function readOpticalData(filename,subdivisions=1) {
    const rows=(await fs.readFile(filename,'utf8')).trim().split(/\r?\n/)
        .map(row=>row.trim().split(/\s+/).map(Number)).filter(row=>row[0]>0);
    const refined=[];
    for(let j=0;j+1<rows.length;j++) {
        const left=rows[j],right=rows[j+1];
        for(let k=0;k<subdivisions;k++) {
            const fraction=k/subdivisions;
            refined.push(left.map((v,c)=>v+(right[c]-v)*fraction));
        }
    }
    refined.push(rows[rows.length-1]);
    return refined.map(([wavenumber,n,k])=>({wavelength:1/wavenumber,n,k})).reverse();
}

const PLANCK_C2=1.438776877;
function planckKernel(wavelength,temperature) {
    const exponent=PLANCK_C2/(wavelength*temperature);
    return exponent>700?0:wavelength**-5/Math.expm1(exponent);
}
export function makeGrainModel(opticalData,radiusUm,solidDensity=3.27,stellarTemperature=3060.752) {
    const count=opticalData.length,wavelengths=opticalData.map(row=>row.wavelength);
    const weights=wavelengths.map((v,j)=>(
        (j+1<count?wavelengths[j+1]-v:0)+(j>0?v-wavelengths[j-1]:0))/2);
    let heating=0,stellarIntegral=0,scatteringIntegral=0,transportIntegral=0,absorptionIntegral=0;
    let minimumAbsorption=Infinity,minimumAsymmetry=Infinity,maximumAsymmetry=-Infinity;
    const coolingWeights=[];
    for(let j=0;j<count;j++) {
        const row=opticalData[j];
        const result=mie([row.n,-row.k],2*Math.PI*radiusUm*1e-4/row.wavelength);
        if(!Number.isFinite(result.absorption) || result.absorption < -1e-9*Math.max(result.extinction,1))
            throw Error('Nonphysical Mie absorption at '+j);
        if(Math.abs(result.asymmetry)>1+1e-10) throw Error('Invalid asymmetry');
        minimumAbsorption=Math.min(minimumAbsorption,result.absorption);
        minimumAsymmetry=Math.min(minimumAsymmetry,result.asymmetry);
        maximumAsymmetry=Math.max(maximumAsymmetry,result.asymmetry);
        const absorption=Math.max(result.absorption,0),stellarWeight=weights[j]*planckKernel(row.wavelength,stellarTemperature);
        stellarIntegral+=stellarWeight;
        scatteringIntegral+=stellarWeight*result.scattering;
        transportIntegral+=stellarWeight*result.scattering*(1-result.asymmetry);
        absorptionIntegral+=stellarWeight*absorption;
        heating+=stellarWeight*absorption;
        coolingWeights.push(weights[j]*absorption);
    }
    const fullIntegral=Math.PI**4/15*(stellarTemperature/PLANCK_C2)**4;
    const coefficient=3/(4*solidDensity*radiusUm*1e-4);
    const cooling=temperature=>{
        let total=0;
        for(let j=0;j<count;j++) total+=coolingWeights[j]*planckKernel(wavelengths[j],temperature);
        return total;
    };
    const dilution=radius=>{const squared=radius**-2;return squared/(2*(1+Math.sqrt(1-squared)));};
    const temperatureAt=radius=>{
        if(radius<1) throw Error('Radius must be at least one stellar radius');
        const target=dilution(radius)*heating;
        let lower=1,upper=stellarTemperature;
        if(cooling(lower)>=target || cooling(upper)<=target) throw Error('Unbracketed grain temperature');
        for(let iteration=0;iteration<40;iteration++) {
            const middle=(lower+upper)/2;
            if(cooling(middle)>target) upper=middle;else lower=middle;
        }
        return (lower+upper)/2;
    };
    const radiusAt=temperature=>{
        const dilutionRequired=cooling(temperature)/heating;
        return dilutionRequired<=0.5 ? 1/(2*Math.sqrt(dilutionRequired*(1-dilutionRequired))) : null;
    };
    return {radiusUm,solidDensity,stellarTemperature,coefficient:coefficient*transportIntegral/fullIntegral,
            rawScatteringCoefficient:coefficient*scatteringIntegral/fullIntegral,
            absorptionCoefficient:coefficient*absorptionIntegral/fullIntegral,
            meanTransportEfficiency:transportIntegral/fullIntegral,
            spectralCoverage:stellarIntegral/fullIntegral,
            minimumAbsorption,minimumAsymmetry,maximumAsymmetry,
            wavelengthRangeCm:[wavelengths[0],wavelengths[count-1]],count,
            temperatureAt,radiusAt,cooling};
}

export async function runAudit(folder,subdivisions=1) {
    const validation=validateMie(),opticalData=await readOpticalData(folder+'/mg2.txt',subdivisions);
    const results=[];
    for(const radius of [0.1,0.3,0.5,1]) {
        const model=makeGrainModel(opticalData,radius);
        results.push({...model,cooling:undefined,temperatureAt:undefined,radiusAt:undefined,
            kappaAtDustRatio1e3:1e-3*model.coefficient,
            dustRatioForPreviousKappa:14.821/model.coefficient,
            temperatures:Object.fromEntries([1,2,2.5,3,4,5,10].map(r=>[r,model.temperatureAt(r)])),
            condensationRadius1100:model.radiusAt(1100),
            condensationRadius1000:model.radiusAt(1000),
            condensationRadius1200:model.radiusAt(1200)});
    }
    const audit={assumptions:{material:'Jena amorphous Mg2SiO4 mg2.txt',solidDensity:3.27,
        stellarTemperature:3060.752,stellarRadiusCm:1.748342799e13,subdivisions,
        illumination:'unattenuated finite-disc blackbody',opacity:'flux-mean Qsca*(1-g), per gram gas'},validation,results};
    await fs.writeFile(folder+'/results_subdivisions_'+subdivisions+'.json',JSON.stringify(audit,null,2));
    return audit;
}

export async function writeTemperatureTable(folder,radiusUm=0.3,subdivisions=4) {
    const data=await readOpticalData(folder+'/mg2.txt',subdivisions);
    const model=makeGrainModel(data,radiusUm);
    const innerTemperature=model.temperatureAt(1),outerTemperature=model.temperatureAt(1024);
    // Uniform log T gives extra spatial resolution close to the photosphere.
    // Inverting the equilibrium relation avoids hundreds of unnecessary root solves.
    const rows=Array.from({length:513},(_,j)=>{
        const temperature=Math.exp(Math.log(innerTemperature)+
            (Math.log(outerTemperature)-Math.log(innerTemperature))*j/512);
        return [j===0?1:j===512?1024:model.radiusAt(temperature),temperature];
    });
    const filename=folder+'/mie_Td_Mg2SiO4_a'+String(radiusUm).replace('.','p')+'um_T3061.dat';
    const header=[
        '# Mg2SiO4 amorphous Jena mg2.txt; grain_radius_um='+radiusUm+'; solid_density_g_cm3=3.27',
        '# Teff_K=3060.752; photospheric_radius_cm=1.748342799e13',
        '# Finite-disc unattenuated blackbody; Qabs-weighted equilibrium; no gas collisions/diffuse radiation',
        '# flux_mean_transport_opacity_cm2_per_g_gas='+model.coefficient*1e-3+' for dust/gas=0.001',
        '# Diagnostic model, not a uniquely calibrated grain population',
        '# r_over_Rstar Td_K'
    ];
    await fs.writeFile(filename,header.join('\n')+'\n'+
        rows.map(row=>row.map(v=>v.toPrecision(15)).join(' ')).join('\n')+'\n');
    return {filename,rows:rows.length,innerTemperature,outerTemperature};
}
