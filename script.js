const arrayOfUnits = [
    ['kg', 'kg', 'mass', 'kilogram', 1],
    ['m', 'm', 'length', 'meter', 1],
    ['s', 's', 'time', 'second', 1],
    ['A', 'A', 'electromag', 'amp', 1],
    ['K', 'K', 'temperature', 'kelvin', 1],
    ['cd', 'cd', 'other', 'candela', 1],
    ['mol', 'mol', 'other', 'mole', 1],
    ['rad', 'rad', 'other', 'radian', 1],
    ['kg', 'mg', 'mass', 'milligram', 1e-6],
    ['kg', 'g', 'mass', 'gram', 1e-3],
    ['kg', 't', 'mass', 'metric ton', 1000],
    ['kg', 'gr', 'mass', 'grain', 64.79891E-6],
    ['kg', 'oz', 'mass', 'ounce', 0.028349523125],
    ['kg', 'lb', 'mass', 'pound', 0.45359237],
    ['kg', 'ton', 'mass', 'short ton', 907.18474],
    ['kg', 'oz.t', 'mass', 'troy ounce', 0.0311034768],
    ['kg', 'Da', 'mass', 'atomic mass unit', 1.6605390666E-27],
    ['m', 'mm', 'length', 'millimeter', 0.001],
    ['m', 'cm', 'length', 'centimeter', 0.01],
    ['m', 'km', 'length', 'kilometer', 1000],
    ['m', 'in', 'length', 'inch', 0.0254],
    ['m', 'ft', 'length', 'feet', 0.3048],
    ['m', 'mi', 'length', 'mile', 1609.344],
    ['m', 'nmi', 'length', 'nautical mile', 1852],
    ['m', 'au', 'length', 'astronomical unit', 149597870700],
    ['m', 'ly', 'length', 'light-year', 9460730472580800E0],
    ['m', 'pc', 'length', 'parsec', 30856775814913673E0],
    ['s', 'ns', 'time', 'nanosecond', 1E-9],
    ['s', '&micro;s', 'time', 'microsecond', 1E-6],
    ['s', 'ms', 'time', 'millisecond', 1E-3],
    ['s', 'min', 'time', 'minute', 60],
    ['s', 'h', 'time', 'hour', 3600],
    ['s', 'd', 'time', 'day', 86400],
    ['s', 'yr', 'time', 'julian year', 31557600],
    ['1/s', 'Hz', 'time', 'hertz', 1],
    ['m2', 'm', 'area', 'sq. meter', 1],
    ['m2', 'ha', 'area', 'hecare', 1E4],
    ['m2', 'mm', 'area', 'sq. millimeter', 1E-6],
    ['m2', 'cm', 'area', 'sq. centimeter', 1E-4],
    ['m2', 'km', 'area', 'sq. kilometer', 1E6],
    ['m2', 'acr', 'area', 'acre', 4046.8564224],
    ['m2', 'in', 'area', 'sq. inch', 0.0254 * 0.0254],
    ['m2', 'ft', 'area', 'sq. feet', 0.3048 * 0.3048],
    ['m2', 'mi', 'area', 'sq. mile', 1609.344 * 1609.344],
    ['m3', 'm', 'volume', 'cubic meter', 1],
    ['m3', 'l', 'volume', 'liter', 1E-3],
    ['m3', 'mm', 'volume', 'cubic millimeter', 1E-9],
    ['m3', 'cm', 'volume', 'cubic centimeter', 1E-6],
    ['m3', 'km', 'volume', 'cubic kilometer', 1E9],
    ['m3', 'fl.oz', 'volume', 'U.S. fluid ounce', 29.5735295625E-9],
    ['m3', 'gal', 'volume', 'U.S. gallon', 3.785411784E-3],
    ['m3', 'in', 'volume', 'cubic inch', 0.0254 * 0.0254 * 0.0254],
    ['m3', 'ft', 'volume', 'cubic feet', 0.3048 * 0.3048 * 0.3048],
    ['m3', 'mi', 'volume', 'cubic mile', 1609.344 * 1609.344 * 1609.344],
    ['m/s', 'm/s', 'speed', 'meter/second', 1],
    ['m/s', 'cm/s', 'speed', 'centimeter/sec', 0.01],
    ['m/s', 'km/h', 'speed', 'kilometer/hour', 1000 / 3600],
    ['m/s', 'km/s', 'speed', 'kilometer/sec', 1000],
    ['m/s', 'ft/s', 'speed', 'feet/second', 0.3048],
    ['m/s', 'mi/h', 'speed', 'mile/hour', 1609.344 / 3600],
    ['m/s', 'mi/s', 'speed', 'mile/second', 1609.344],
    ['m/s', 'c', 'speed', 'light speed', 299792458],
    ['kg m/s2', 'N', 'force', 'newton', 1],
    ['kg m/s2', 'dyn', 'force', 'dyne', 1E-5],
    ['kg m/s2', 'lbf', 'force', 'pound force', 4.448222],
    ['kg/m s2', 'Pa', 'pressure', 'pascal', 1],
    ['kg/m s2', 'kPa', 'pressure', 'kilopascal', 1000],
    ['kg/m s2', 'bar', 'pressure', 'bar', 1E5],
    ['kg/m s2', 'mbar', 'pressure', 'millibar', 100],
    ['kg/m s2', 'psi', 'pressure', 'lfb/in&super2;', 6894.757],
    ['K', '°C', 'temperature', 'deg. Celcuis', 1, -273.15],
    ['K', '°F', 'temperature', 'deg. Fahrenheit', 5 / 9, -459.67],
    ['kg m2/s2', 'J', 'energy', 'joule', 1],
    ['kg m2/s2', 'kW.h', 'energy', 'kilowatt-hour', 3.6E6],
    ['kg m2/s2', 'BTU', 'energy', 'ISO BTU', 1055.6],
    ['kg m2/s2', 'erg', 'energy', 'joule', 1E-7],
    ['kg m2/s2', 'kcal', 'energy', 'kilocalorie', 4184],
    ['kg m2/s2', 'tTNT', 'energy', 'ton-TNT', 4.184E9],
    ['kg m2/s2', 'ft.lbf', 'energy', 'foot pound force', 1.355818],
    ['kg m2/s2', 'eV', 'energy', 'electronvolt', 1.602176634E-19],
    ['kg m2/s3', 'W', 'power', 'watt', 1],
    ['kg m2/s3', 'kW', 'power', 'kilowatt', 1E3],
    ['kg m2/s3', 'MW', 'power', 'megawatt', 1E6],
    ['kg m2/s3', 'hp', 'power', 'horsepower', 745.69987158227],
    ['s A', 'C', 'electromag', 'coulomb', 1],
    ['kg m2/s3 A', 'V', 'electromag', 'volt', 1],
    ['s4 A2/kg m2', 'F', 'electromag', 'farad', 1],
    ['kg m2/s3 A2', 'Ω', 'electromag', 'ohm', 1],
    ['s3 A2/kg m2', 'S', 'electromag', 'siemens', 1],
    ['kg m2/s2 A', 'Wb', 'electromag', 'weber', 1],
    ['kg/s2 A', 'T', 'electromag', 'tesla', 1],
    ['kg m2/s2 A2', 'H', 'electromag', 'henry', 1],
    ['cd/m2', 'lx', 'other', 'lux', 1],
    ['rad', 'deg', 'other', 'degree', Math.PI / 180],
];
const arrayOfConstants = [
    ['1', Math.PI, 'π', 'basic', 'pi', 'pi'],
    ['m/s', 299792458, 'c', 'electromag', 'speed of light'],
    ['m3/kg s2', 6.6743E-11, 'G', 'astronomy', 'newtonian gravitational constant'],
    ['kg m2/s', 6.62607015E-34, 'h', 'particle', 'planck constant'],
    ['kg m2/s', 1.054571817E-34, 'ħ', 'particle', 'reduced planck constant', 'hbar'],
    ['A s', 1.602176634E-19, 'e', 'particle', 'elementry charge'],
    ['m/s2 A2', 1.25663706212E-6, 'µ₀', 'electromag', 'vacuum magnetic permittivity', 'mu0'],
    ['s4 A2/m3 kg', 8.8541878128E-12, 'ε₀', 'electromag', 'vacuum electric permittivity', 'epsilon0'],
    ['s2 A/kg m2', 483597.9E9, 'Kj', 'electromag', 'Josephson constant'],
    ['kg m2/s3 A2', 25812.807, 'Rk', 'electromag', 'von Klitzing constant'],
    ['kg m2/s2 A', 2.067833848E-15, 'Φ₀', 'electromag', 'magnetic flux quantum', 'phi0'],
    ['s3 A2/kg m2', 7.748091729E-5, 'G₀', 'electromag', 'conductance quantum', 'G0'],
    ['kg', 9.1093837015E-31, 'Me', 'particle', 'electron mass'],
    ['kg', 1.67262192369E-27, 'Mp', 'particle', 'proton mass'],
    ['1', 1836.15267343, 'Mp/e', 'particle', 'proton to electron mass ratio', 'Mp2Me'],
    ['1', 7.2973525693E-3, 'α', 'particle', 'fine structure constant', 'alpha'],
    ['1', 137.035999084, '1/α', 'particle', 'inverse fine structure constant', 'alpha1'],
    ['1/s', 3.2898419602508E15, 'cR∞', 'particle', 'Rydberg frequency', 'cRinfity'],
    ['kg m2/s2 K', 1.380649E-23, 'k', 'thermodyn', 'Boltzmann constant'],
    ['1/mol', 6.02214076E23, 'L', 'particle', 'Avogadro constant'],
    ['kg m2/s2 K mol', 8.31446261815324, 'R', 'thermodyn', 'molar gas constant'],
    ['s A/mol', 96485.33212, 'F', 'electromag', 'Faraday constant'],
    ['kg/s3 K4', 5.670374419E-8, 'σ', 'thermodyn', 'Stefan-Boltzmann constant', 'sigma'],
    ['kg m2/s2', 1.602176634E-19, 'eV', 'electromag', 'electronvolt'],
    ['kg', 1.6605390666E-27, 'u', 'particle', 'unified mass unit'],
    ['kg', 1.67492749804E-27, 'Mn', 'particle', 'neutron mass'],
    ['kg m3/s4 A2', 8.9875517923E9, 'ke', 'electromag', 'Coulomb constant'],
    ['kg m2/s3 A2', 376.730313668, 'Z₀', 'electromag', 'characteristic impedance of vacuum', 'Z0'],
    ['m2/s2', 333550, '♆HF', 'thermodyn', 'heat of fusion by mass of water', 'h2oHF'],
    ['m2/s2', 2257000, '♆HV', 'thermodyn', 'heat of vaporization by mass of water', 'h2oHV'],
    ['m2/s2 K', 4184, '♆HC', 'thermodyn', 'heat capacity of water by mass at 20C', 'h2oHC'],
    ['m/s2', 9.80665, 'g₀', 'astronomy', 'standard gravity', 'g0'],
    ['kg', 1.98847E30, 'M⊙', 'astronomy', 'mass of Sun', 'Msun'],
    ['kg', 5.9722E24, 'M🜨', 'astronomy', 'mass of Earth', 'Mearth'],
    ['m', 149597870700, 'au', 'astronomy', 'astronomical unit length'],
    ['m', 6378137, 'R🜨', 'astronomy', 'equatorial radius of Earth', 'Rearth'],
];
const arrayOfFormula = [
    ['ohms law', 'electrical',
        ['v = r * i', 'i = v / r', 'r = v / i'],
        ['v=kg m2/s3 A', 'i=A', 'r=kg m2/s3 A2']],
    ['watts law', 'electrical',
        ['p = v * i', 'i = p / v', 'v = p / i'],
        ['p=kg m2/s3', 'i=A', 'v=kg m2/s3 A']],
    ['ohms/watts law', 'electrical',
        ['p = r * i**2', 'r = p / i**2', 'i = sqrt(p / i)'],
        ['p=kg m2/s3', 'r=kg m2/s3 A2', 'i=A']],
    ['ohms/watts law', 'electrical',
        ['r = v**2 / p', 'v = sqrt(r * p)', 'p = v**2 / r'],
        ['r=kg m2/s3 A2', 'v=kg m2/s3 A', 'p=kg m2/s3']],
    ['sphere volume, radius', 'geometry',
        ['v = 4/3*PI*r**3', 'r = cbrt(v*3 / (4*PI))'],
        ['v=m3', 'r=m']],
    ['sphere area, radius ', 'geometry',
        ['a = 4*PI * r**2', 'r = sqrt(a*PI / 4)'],
        ['a=m2', 'r=m']],
    ['sphere volume, diameter', 'geometry',
        ['v = PI/6 * d**3', 'd = cbrt(6*v / PI)'],
        ['v=m3', 'd=m']],
    ['sphere area, diameter', 'geometry',
        ['a = PI * d**2', 'd = sqrt(a/PI)'],
        ['a=m2', 'd=m']],
    ['newton 2nd law', 'mechanics',
        ['F = m*a', 'm = F / a', 'a = F / m'],
        ['F=kg m/s2', 'm=kg', 'a=m/s2']],
    ['accelerated motion*', 'mechanics',
        ['a = (v-v0) / delta_t', 'delta_t = (v-v0) / a', 'v = a*delta_t + v0', 'v0 = v - a*delta_t'],
        ['a=m/s2', 'v=m/s', 'v0=m/s', 'delta_t=s']],
    ['accelerated motion', 'mechanics',
        ['x = 1/2*a*t**2 + v*t', 't = (sqrt(2*a*x + v**2) - v) / a', 't = -(sqrt(2*a*x + v**2) + v) / a', 'v = x/t - a*t/2', 'a = (2*x - 2*t*v)/t**2'],
        ['x=m', 'a=m/s2', 'v=m/s', 't=s']],
    ['circular motion', 'mechanics',
        ['a = v**2 / r', 'r = v**2 / a', 'v = sqrt(a*r)'],
        ['a=m/s2', 'v=m/s', 'r=m']],
    ['circular motion', 'mechanics',
        ['F = m*v**2/r', 'r = m*v**2 / F', 'v = sqrt(F*r/m)', 'm = F*r / v**2'],
        ['F=kg m/s2', 'm=kg', 'v=m/s', 'r=m']],
    ['gravitational potential energy', 'mechanics',
        ['U = -$G*M*m / r', 'r = -$G*M*m / U', 'm = -U*r / $G*M', 'M = -U*r / $G*n'],
        ['U=kg m2/s2', 'M=kg', 'm=kg', 'r=m']],
    ['potential energy, lifted mass', 'mechanics',
        ['delta_U = m*$g0*h', 'm = delta_U / ($g0*h)', 'h = delta_U / (m*$g0)'],
        ['delta_U=kg m2/s2', 'm=kg', 'h=m']],
    ['elastic potential energy', 'mechanics',
        ['delta_u = (1/2) * k * x**2', 'k = 2*delta_u / x**2', 'x = sqrt(2*delta_u / k)', 'x = -sqrt(2*delta_u / k)'],
        ['delta_u=kg m2/s2', 'k=1', 'x=m']],
    ['kinetic energy', 'mechanics',
        ['E = (1/2) * m * v**2', 'm = 2*E / v**2', 'v = sqrt(2*E / m)', 'v = -sqrt(2*E / m)'],
        ['E=kg m2/s2', 'm=kg', 'v=m/s']],
    ['gravitational attraction', 'mechanics',
        ['F = $G*m1*m2/r**2', 'r = sqrt($G*m1*m2 / F)', 'm1 = F*r**2 / ($G*m2)', 'm2 = F*r**2 / ($G*m1)'],
        ['F=kg m/s2', 'm1=kg', 'm2=kg', 'r=m']],
    ['gravitational attraction', 'mechanics',
        ['a = $G*m/r**2', 'r = sqrt($G*m / a)', 'm = a*r**2 / $G'],
        ['a=m/s2', 'm=kg', 'r=m']],
    ['rocket equation*', 'mechanics',
        ['delta_v = ve * log(m0 / m1)', 've = log(m0 / m1) /  delta_v', 'm0 = m1 * exp(v / delta_v)', 'm1 = m0 * exp(-v / delta_v)'],
        ['delta_v=m/s', 've=m/s', 'm0=kg', 'm1=kg']],
    ['gas law*', 'thermodynamics',
        ['p1 = p2*t1*v2 / t2*v1', 'v1 = v2*p2*t1 / p1*t2', 't1 = t2*p1*v1 / p2*v2'],
        ['t1=K', 'p1=kg/m s2', 'v1=m3', 't2=K', 'p2=kg/m s2', 'v2=m3']],
    ['gas law, constant volume*', 'thermodynamics',
        ['p1 = p2*t1 / t2', 't1 = t2*p1 / p2'],
        ['t1=K', 'p1=kg/m s2', 't2=K', 'p2=kg/m s2']],
    ['gas law, constant pressure*', 'thermodynamics',
        ['v1 = v2*t1 / t2', 't1 = t2*v1 / v2'],
        ['t1=K', 'v1=m3', 't2=K', 'v2=m3']],
    ['gas law, constant temperature*', 'thermodynamics',
        ['p1 = p2*v2 / v1', 'v1 = v2*p2 / p1'],
        ['p1=kg/m s2', 'v1=m3', 'p2=kg/m s2', 'v2=m3']],
    ['ideal gas law', 'thermodynamics',
        ['p = n*$R*t / v', 'v = n*$R / p', 'n = p*v / ($R*t)', 't = p*v / (n*$R)'],
        ['t=K', 'p=kg/m s2', 'n=mol', 'v=m3']],
];
const arrayOfButtons = [
    [
        ['CE/C', 'CE/C', 'CE/C', 'CE/C'],
        ['π', 'π', 'm', 'π'],
        ['LN', 'LOG', 'mm', 'c'],
        ['EE', 'EE', 'km', 'µ₀'],
        ['(', '(', 'ft', 'ε₀'],
        [')', ')', 'in', 'G₀'],
        ['÷', '÷', 'mi', 'Φ₀'],
        ['¹/ₓ', '¹/ₓ', '1/UN', 'F']
    ], [
        ['UNIT', 'UNIT', 'UNIT', ''],
        ['STR', 'SUM', 'kg', 'k'],
        ['eˣ', '10ˣ', 'g', 'R'],
        ['7', '7', 'mg', 'σ'],
        ['8', '8', 'lb', '♆HF'],
        ['9', '9', 'oz', '♆HV'],
        ['×', '×', 'l', '♆HC'],
        ['yˣ', 'ˣ√y', 'gal', 'eV']
    ], [
        ['CNST', 'CNST', '', 'CNST'],
        ['RCL', 'EXCH', 's', 'h'],
        ['SIN', 'SIN⁻¹', 'min', 'ħ'],
        ['4', '4', 'h', 'e'],
        ['5', '5', 'd', 'cR∞'],
        ['6', '6', 'yr', 'L'],
        ['-', '-', 'rad', 'u'],
        ['√x', '³√x', 'deg', 'ke']
    ], [
        ['2nd', '2nd', 'KWN', ''],
        ['HYP', 'HYP', 'A', 'Me'],
        ['COS', 'COS⁻¹', 'V', 'Mp'],
        ['1', '1', 'Ω', 'Mp/e'],
        ['2', '2', 'J', 'Mn'],
        ['3', '3', 'W', 'α'],
        ['+', '+', 'N', '1/α'],
        ['x²', 'x³', 'Pa', 'Z₀']
    ], [
        ['FLT+', 'FLT+', 'UNKN', ''],
        ['RAD+', 'RAD+', 'K', 'G'],
        ['TAN', 'TAN⁻¹', '°C', 'g₀'],
        ['0', '0', '°F', 'M⊙'],
        ['.', '.', 'cd', 'M🜨'],
        ['±', '±', 'mol', 'au'],
        ['=', '=', '#', 'R🜨'],
        ['LIST', 'LIST', 'LIST', 'LIST']
    ]
];
var buttonElements = new Map([
    ['UNIT', undefined],
    ['LIST', undefined],
    ['CNST', undefined],
    ['2nd', undefined],
    ['1/UN', undefined],
    ['HYP', undefined],
    ['RAD+', undefined],
    ['FLT+', undefined],
]);
const buttonColorKeys = [
    ['CE/C', 'MODE', 'CNST', 'UNIT', '2nd', 'LIST', 'KWN', 'UNKN', '1/UN', '#', 'RAD+', 'FLT+', 'HYP'],
    ['.', '±', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'EE'],
    ['+', '-', '×', '÷', '=', '(', ')'],
    ['π', 'STR', 'RCL', 'SIN', 'COS', 'TAN', 'LN', 'eˣ', 'yˣ', '¹/ₓ', '√x', 'x²', 'LOG',
        ,
        'SUM', 'EXCH', '10ˣ', 'ˣ√y', '³√x', 'x³', 'SIN⁻¹', 'COS⁻¹', 'TAN⁻¹']
];
function hsl(h, s, l) {
    var a = s * Math.min(l, 1 - l);
    function f(n) {
        var k = (n + h / 30) % 12;
        var color = l - a * Math.max(Math.min(k - 3, 9 - k, 1), -1);
        return ('0' + Math.round(255 * color).toString(16)).slice(-2);
    }
    ;
    return `#${f(0)}${f(8)}${f(4)}`;
}
const buttonColors = [hsl(30, 1, .2), hsl(210, 1, .2), hsl(90, 1, .15), hsl(255, 1, .15), hsl(165, 1, .2), hsl(330, 1, .2), hsl(0, 0, .5)];
function format(n) {
    if (!isFinite(n))
        return n.toString();
    const significantDigits = 8;
    var components = ((n.toExponential()).toString()).split('e');
    var mantisa = parseFloat(components[0]);
    mantisa = parseFloat(mantisa.toFixed(significantDigits - 1));
    var exponent = parseInt(components[1]);
    if (fltSciEng == 'FLT+' && exponent < significantDigits && exponent >= -3) {
        var digits = significantDigits - 1;
        if (exponent < 0)
            digits -= exponent;
        return parseFloat(n.toFixed(digits)).toString();
    }
    else {
        if (fltSciEng == 'ENG+') {
            var fix = exponent % 3;
            if (fix == -1) {
                mantisa *= 100;
                exponent -= 2;
            }
            else if (fix == -2) {
                mantisa *= 10;
                exponent -= 1;
            }
            else {
                mantisa *= 10 ** fix;
                exponent -= fix;
            }
        }
        return mantisa.toString() + '×10<sup><small>' + exponent.toString() + '</small></sup>';
    }
}
class ShortUnits {
    constructor(str, powers) {
        if (str == undefined) {
            this.powers = powers.slice();
            return;
        }
        this.powers = [0, 0, 0, 0, 0, 0, 0, 0];
        var fraction = str.split('/');
        for (var i = 0; i < fraction.length; i++) {
            if (fraction[i] == '1')
                continue;
            var multiple = fraction[i].split(' ');
            var j;
            for (j = 0; j < multiple.length; j++) {
                var name = multiple[j].match(/[a-zA-Z]+/)[0];
                var inx = 0;
                for (inx = 0; inx < ShortUnits.names.length; inx++)
                    if (name == ShortUnits.names[inx])
                        break;
                var exps = multiple[j].match(/[0-9]+/);
                var exp = 1;
                if (exps != null && exps.length > 0 && exps[0] != '')
                    exp = parseFloat(exps[0]);
                if (i > 0)
                    exp = -exp;
                this.powers[inx] = exp;
            }
        }
    }
    copy() {
        return new ShortUnits(undefined, this.powers);
    }
    isEqual(other) {
        return (this.powers.toString() == other.powers.toString());
    }
    isEqualArray(other) {
        return (this.powers.toString() == other.toString());
    }
    toString() {
        var ret = '';
        for (var i = 0; i < 8; i++) {
            if (this.powers[i] > 1)
                ret += ShortUnits.names[i] + this.powers[i].toString() + ' ';
            else if (this.powers[i] == 1)
                ret += ShortUnits.names[i] + ' ';
        }
        ret = ret.trim();
        if (ret == '')
            ret = '1';
        ret += '/';
        for (var i = 0; i < 8; i++) {
            if (this.powers[i] < -1)
                ret += ShortUnits.names[i] + (-this.powers[i].toString()) + ' ';
            else if (this.powers[i] == -1)
                ret += ShortUnits.names[i] + ' ';
        }
        ret = ret.trim();
        if (ret == '1/')
            ret = '';
        if (ret.substr(-1, 1) == '/')
            ret = ret.substr(0, ret.length - 1);
        return ret;
    }
}
ShortUnits.names = ['kg', 'm', 's', 'K', 'A', 'cd', 'mol', 'rad'];
class Unit {
    constructor(shortUnits, name, group, desc, factor, offset = 0) {
        this.units = new ShortUnits(shortUnits);
        this.name = name;
        this.group = group;
        this.desc = desc;
        this.factor = factor;
        this.offset = offset;
    }
    nPowers() {
        var n = 0;
        for (var i = 0; i < 8; i++)
            if (this.units.powers[i] != 0)
                n++;
        return n;
    }
    isComposite() {
        return this.name.includes('/');
    }
    isComplex() {
        return ((this.nPowers() > 1 || this.units.powers[this.index()] != 1 || this.name.includes('°')) && !this.isComposite());
    }
    index() {
        for (var i = 0; i < 8; i++)
            if (this.units.powers[i] != 0)
                return i;
        return -1;
    }
    names() {
        var un = ['', '', '', '', '', '', '', ''];
        if (this.isComplex()) {
            for (var i = 0; i < 8; i++)
                if (this.units.powers[i] != 0)
                    un[i] = ShortUnits.names[i];
        }
        else {
            var terms = this.name.split('/');
            for (var i = 0; i < terms.length; i++)
                un[findUnitByName(terms[i]).index()] = terms[i];
        }
        return un;
    }
}
var listOfUnits = [];
function findUnitByName(name) {
    for (var i = 0; i < listOfUnits.length; i++)
        if (name == listOfUnits[i].name)
            return listOfUnits[i];
    return undefined;
}
function findUnit(powers, name) {
    for (var i = 0; i < listOfUnits.length; i++)
        if (listOfUnits[i].units.isEqualArray(powers) && (name == undefined || name == listOfUnits[i].name))
            return listOfUnits[i];
    return undefined;
}
class Constant {
    constructor(shortUnits, value, name, group, desc, altName = undefined) {
        this.units = new ShortUnits(shortUnits);
        this.value = value;
        this.name = name;
        this.group = group;
        this.desc = desc;
        if (altName == undefined)
            this.altName = name;
        else
            this.altName = altName;
    }
    toMeasure() {
        return new Measurement(this.value, this.units.powers);
    }
}
var listOfConstants = [];
class Measurement {
    constructor(value, unitPowers, unitNames, complex, formulaVar) {
        this.value = value;
        this.complexUnits = complex;
        this.formulaVar = formulaVar;
        if (unitPowers == undefined)
            this.unitPowers = [0, 0, 0, 0, 0, 0, 0, 0];
        else
            this.unitPowers = unitPowers.slice();
        if (unitNames == undefined) {
            this.unitNames = [];
            for (var i = 0; i < 8; i++)
                if (this.unitPowers[i] == 0)
                    this.unitNames.push('');
                else
                    this.unitNames.push(ShortUnits.names[i]);
        }
        else
            this.unitNames = unitNames.slice();
    }
    copy() {
        return new Measurement(this.value, this.unitPowers, this.unitNames, this.complexUnits, this.formulaVar);
    }
    factor() {
        var ret = 1;
        var unit = findUnit(this.unitPowers, this.complexUnits);
        if (unit)
            return unit.factor;
        for (var i = 0; i < 8; i++) {
            if (this.unitPowers[i] != 0) {
                var searchUnit = [0, 0, 0, 0, 0, 0, 0, 0];
                searchUnit[i] = 1;
                unit = findUnit(searchUnit, this.unitNames[i]);
                if (unit)
                    ret *= unit.factor ** this.unitPowers[i];
            }
        }
        return ret;
    }
    toDisplayForm() {
        const superscripts = ['', '', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
        var unit;
        if (this.isScalar() && currentMode != 'mode-unit')
            return { value: format(this.value), numerator: '', denominator: '' };
        unit = findUnit(this.unitPowers, this.complexUnits);
        if (unit)
            return { value: format(this.value / unit.factor + unit.offset), numerator: unit.name, denominator: '&nbsp;' };
        var factor = 1;
        var numerator = '';
        var denominator = '';
        for (var i = 0; i < 8; i++) {
            if (this.unitPowers[i] != 0) {
                var searchUnit = [0, 0, 0, 0, 0, 0, 0, 0];
                searchUnit[i] = 1;
                unit = findUnit(searchUnit, this.unitNames[i]);
                if (unit) {
                    if (this.unitPowers[i] > 0)
                        numerator += unit.name + superscripts[this.unitPowers[i]] + ' ';
                    else
                        denominator += unit.name + superscripts[-this.unitPowers[i]] + ' ';
                    factor *= unit.factor ** this.unitPowers[i];
                }
            }
        }
        if (currentMode == 'mode-unit') {
            if (unitSign == 1)
                numerator += '◆';
            else
                denominator += '◆';
        }
        if (denominator != '' && numerator == '')
            numerator = '1';
        numerator = numerator.trim();
        denominator = denominator.trim();
        if (numerator != '' && denominator == '')
            denominator = '&nbsp;';
        return { value: format(this.value / factor), numerator: numerator, denominator: denominator };
    }
    nPowers() {
        var n = 0;
        for (var i = 0; i < 8; i++)
            if (this.unitPowers[i] != 0)
                n++;
        return n;
    }
    isScalar() {
        return this.unitPowers.toString() == [0, 0, 0, 0, 0, 0, 0, 0].toString();
    }
    isScalarOrRadian() {
        var pows = this.unitPowers.slice();
        pows[7] = 0;
        return pows.toString() == [0, 0, 0, 0, 0, 0, 0, 0].toString();
    }
}
class Formulas {
    constructor(description, group, solutions, varShortUnits) {
        this.desc = description;
        this.group = group;
        this.solutions = solutions;
        this.varNames = [];
        this.varUnits = [];
        if (varShortUnits)
            for (var i = 0; i < varShortUnits.length; i++) {
                var v = varShortUnits[i].split('=');
                this.varNames[i] = v[0];
                this.varUnits[i] = new ShortUnits(v[1]);
            }
    }
    copy() {
        var f = new Formulas(this.desc, this.group, this.solutions.slice());
        f.matching = this.matching;
        f.varNames = this.varNames.slice();
        for (var i = 0; i < this.varUnits.length; i++)
            f.varUnits[i] = this.varUnits[i].copy();
        return f;
    }
    matchVariables() {
        for (var i = 0; i < knowns.length; i++)
            knowns[i].formulaVar = '';
        for (var i = 0; i < this.varUnits.length; i++) {
            for (var j = 0; j < knowns.length; j++)
                if (this.varUnits[i].isEqualArray(knowns[j].unitPowers)) {
                    if (knowns[j].formulaVar == '') {
                        knowns[j].formulaVar = this.varNames[i];
                        break;
                    }
                }
        }
    }
    findMatchsForUnknown() {
        var solutions = [];
        for (var i = 0; i < this.solutions.length; i++) {
            var formula = this.solutions[i];
            var formulaVar = formula.split('=')[0].trim();
            if (formulaVar == knowns[0].formulaVar)
                solutions.push(formula);
        }
        return solutions;
    }
    solve() {
        var java = this.matching;
        java = java.replace(/PI/g, 'Math.PI');
        java = java.replace(/exp/g, 'Math.exp');
        java = java.replace(/log/g, 'Math.log');
        java = java.replace(/sin/g, 'Math.sin');
        java = java.replace(/cos/g, 'Math.cos');
        java = java.replace(/tan/g, 'Math.tan');
        java = java.replace(/asin/g, 'Math.asin');
        java = java.replace(/acos/g, 'Math.acos');
        java = java.replace(/atan/g, 'Math.atan');
        java = java.replace(/sqrt/g, 'Math.sqrt');
        java = java.replace(/cbrt/g, 'Math.cbrt');
        var str = "'use strict'; ";
        for (var i = 0; i < listOfConstants.length; i++)
            str += 'var $' + listOfConstants[i].altName + ' = ' + listOfConstants[i].value + '; ';
        str += 'function root4(x) { return x**(1/4) } ';
        for (var i = 1; i < knowns.length; i++)
            if (knowns[i].formulaVar != '')
                str += 'var ' + knowns[i].formulaVar + ' = ' + knowns[i].value + '; ';
        str += 'var ' + java + '; ';
        str += 'return ' + knowns[0].formulaVar;
        return Function(str)();
    }
    prettyMatching() {
        var pretty = this.matching;
        for (var i = 0; i < listOfConstants.length; i++)
            if (listOfConstants[i].name != listOfConstants[i].altName) {
                var regex = new RegExp('$' + listOfConstants[i].altName, 'g');
                pretty = pretty.replace(regex, listOfConstants[i].name);
            }
        pretty = pretty.replace(/$/g, '');
        pretty = pretty.replace(/PI/g, 'π');
        pretty = pretty.replace(/asin/g, 'sin⁻¹');
        pretty = pretty.replace(/acos/g, 'cos⁻¹');
        pretty = pretty.replace(/atan/g, 'tan⁻¹');
        pretty = pretty.replace(/\*\*2/g, '²');
        pretty = pretty.replace(/\*\*3/g, '³');
        pretty = pretty.replace(/\*\*4/g, '⁴');
        pretty = pretty.replace(/\*\*5/g, '⁵');
        pretty = pretty.replace(/sqrt/g, '√');
        pretty = pretty.replace(/cbrt/g, '∛');
        pretty = pretty.replace(/cbrt/g, '∛');
        pretty = pretty.replace(/root4/g, '∜');
        pretty = pretty.replace(/delta_/g, 'Δ');
        return pretty;
    }
}
var listOfFormulas = [];
function findMatchingFormulas(type) {
    var sortGrouping;
    var found = [];
    function swapFormulaVars(f, i, j) {
        var n = f.varNames[i];
        f.varNames[i] = f.varNames[j];
        f.varNames[j] = n;
        var t = f.varUnits[i].copy();
        f.varUnits[i] = f.varUnits[j].copy();
        f.varUnits[j] = t;
        var g = sortGrouping[i];
        sortGrouping[i] = sortGrouping[j];
        sortGrouping[j] = g;
    }
    function permutateVarsAndAddToFound(formula, first, last) {
        if (last == first + 1)
            found.push(formula.copy());
        else {
            permutateVarsAndAddToFound(formula.copy(), first, last - 1);
            for (var i = first; i < last - 1; i++) {
                if (((last - first) & 1) == 0) {
                    if (sortGrouping[i] == sortGrouping[last - 1]) {
                        swapFormulaVars(formula, i, last - 1);
                        permutateVarsAndAddToFound(formula.copy(), first, last - 1);
                    }
                }
                else {
                    if (sortGrouping[first] == sortGrouping[last - 1]) {
                        swapFormulaVars(formula, first, last - 1);
                        permutateVarsAndAddToFound(formula.copy(), first, last - 1);
                    }
                }
            }
        }
    }
    function sortRepeatedVarsToEnd(formula, exclude) {
        sortGrouping = (new Array(formula.varUnits.length)).fill('');
        for (var i = 0; i < formula.varUnits.length; i++)
            for (var j = i + 1; j < formula.varUnits.length; j++)
                if (formula.varUnits[i].isEqual(formula.varUnits[j]) && formula.varNames[i] != exclude && formula.varNames[j] != exclude) {
                    sortGrouping[i] = formula.varUnits[i].toString();
                    sortGrouping[j] = formula.varUnits[i].toString();
                }
        var n = formula.varUnits.length;
        var didSwap;
        do {
            didSwap = false;
            for (i = 1; i < n; i++) {
                if (sortGrouping[i - 1] > sortGrouping[i]) {
                    swapFormulaVars(formula, i - 1, i);
                    didSwap = true;
                }
            }
            n--;
        } while (didSwap);
        for (i = 0; i < formula.varUnits.length; i++)
            if (sortGrouping[i] != '')
                break;
        return i;
    }
    for (var i = 0; i < listOfFormulas.length; i++) {
        var formula = listOfFormulas[i];
        var knownMatches = (new Array(knowns.length)).fill(false);
        for (var j = 0; j < formula.varUnits.length; j++) {
            for (var k = 0; k < knowns.length; k++)
                if (formula.varUnits[j].isEqualArray(knowns[k].unitPowers) && knownMatches[k] == false) {
                    knownMatches[k] = true;
                    break;
                }
        }
        var formulaMatches = (new Array(formula.varUnits.length)).fill(false);
        for (var j = 0; j < formula.varUnits.length; j++) {
            for (var k = 0; k < knowns.length; k++)
                if (formula.varUnits[j].isEqualArray(knowns[k].unitPowers) && formulaMatches[j] == false) {
                    formulaMatches[j] = true;
                    break;
                }
        }
        var match = false;
        switch (type) {
            case 'missing':
                match = (!knownMatches.includes(false) && formulaMatches.includes(false));
                break;
            case 'extra':
                match = (knownMatches.includes(false) && !formulaMatches.includes(false));
                break;
            case 'exact':
                match = (!knownMatches.includes(false) && !formulaMatches.includes(false));
                break;
        }
        if (match) {
            formula.matchVariables();
            var matchSubFormulas = formula.findMatchsForUnknown();
            for (var j = 0; j < matchSubFormulas.length; j++) {
                var f = formula.copy();
                f.matching = matchSubFormulas[j];
                var exclude = (f.matching.split('='))[0].trim();
                if (/\*$/.test(f.desc))
                    permutateVarsAndAddToFound(f, sortRepeatedVarsToEnd(f, exclude), f.varNames.length);
                else
                    found.push(f);
            }
        }
    }
    return found;
}
var mantisa = '';
var exponent = '';
var entryMode = 'number';
var operators = [];
var operands = [];
var memory = new Measurement(0);
var unitSign = 1;
var currentMode;
var hyper = false;
var degRad = 'RAD+';
var fltSciEng = 'FLT+';
var exactFormulas;
var implicitFormula;
var missingTermFormulas;
var extraTermFormulas;
var listedFormulas = implicitFormula;
var inxFormulas = 0;
var knowns = [new Measurement(0)];
window.onload = function () {
    function create(constructor, argList) {
        return new (Function.prototype.bind.apply(constructor, [null].concat(argList)));
    }
    listOfConstants = [];
    for (var i = 0; i < arrayOfConstants.length; i++)
        listOfConstants.push(create(Constant, arrayOfConstants[i]));
    listOfUnits = [];
    for (var i = 0; i < arrayOfUnits.length; i++)
        listOfUnits.push(create(Unit, arrayOfUnits[i]));
    listOfFormulas = [];
    for (var i = 0; i < arrayOfFormula.length; i++)
        listOfFormulas.push(create(Formulas, arrayOfFormula[i]));
    fillTreeInHTML(listOfConstants, 'constTreeView', 'selectConstByName');
    fillTreeInHTML(listOfUnits, 'unitTreeView', 'selectUnitByName');
    clearButton(true);
    wireUpTreeTogglerInHTML();
    setButtonMode('mode-norm');
    setupScroll();
};
var landscapeMediaQuery = window.matchMedia('only screen and (orientation: landscape)');
landscapeMediaQuery.addEventListener('change', function () { setButtonMode(currentMode); });
function setButtonMode(newMode) {
    var modes = ['mode-norm', 'mode-2nd', 'mode-unit', 'mode-const'];
    var mInx = modes.indexOf(newMode);
    for (var [key] of buttonElements)
        buttonElements.set(key, undefined);
    var groups = [];
    for (k = 0; k < listOfConstants.length; k++)
        if (!(groups.includes(listOfConstants[k].group)))
            groups.push(listOfConstants[k].group);
    for (var i = 0; i < 5; i++)
        for (var j = 0; j < 8; j++) {
            var button;
            if (landscapeMediaQuery.matches)
                button = document.getElementById('x' + j.toString() + 'y' + i.toString());
            else {
                if ((newMode == 'mode-norm' || newMode == 'mode-2nd') && j >= 3)
                    button = document.getElementById('x' + (j - 3).toString() + 'y' + (i + 3).toString());
                else
                    button = document.getElementById('x' + i.toString() + 'y' + j.toString());
            }
            button.innerHTML = arrayOfButtons[i][j][mInx];
            if (newMode == 'mode-norm' || newMode == 'mode-2nd') {
                for (var k = 0; k < buttonColorKeys.length; k++)
                    if (buttonColorKeys[k].includes(button.innerHTML)) {
                        button.style.backgroundColor = buttonColors[k];
                        break;
                    }
                if (k == buttonColorKeys.length)
                    button.style.backgroundColor = buttonColors[6];
            }
            else if (newMode == 'mode-unit') {
                var unit = findUnitByName(button.innerHTML);
                if (buttonColorKeys[0].includes(button.innerHTML))
                    button.style.backgroundColor = buttonColors[0];
                else if (unit) {
                    if (unit.units.isEqual(new ShortUnits('m')))
                        button.style.backgroundColor = buttonColors[1];
                    else if (unit.units.isEqual(new ShortUnits('kg')))
                        button.style.backgroundColor = buttonColors[2];
                    else if (unit.units.isEqual(new ShortUnits('s')))
                        button.style.backgroundColor = buttonColors[3];
                    else if (unit.units.isEqual(new ShortUnits('K')) || unit.units.isEqual(new ShortUnits('A')) ||
                        unit.units.isEqual(new ShortUnits('mol')) || unit.units.isEqual(new ShortUnits('cd')) ||
                        unit.units.isEqual(new ShortUnits('rad')))
                        button.style.backgroundColor = buttonColors[4];
                    else
                        button.style.backgroundColor = buttonColors[5];
                }
                else
                    button.style.backgroundColor = buttonColors[6];
            }
            else if (newMode == 'mode-const') {
                var cnst = findConstantByName(button.innerHTML);
                if (buttonColorKeys[0].includes(button.innerHTML))
                    button.style.backgroundColor = buttonColors[0];
                else if (cnst) {
                    var inx = groups.indexOf(cnst.group);
                    if (inx >= 4)
                        button.style.backgroundColor = buttonColors[5];
                    else if (inx < 0)
                        button.style.backgroundColor = buttonColors[6];
                    else
                        button.style.backgroundColor = buttonColors[inx + 1];
                }
                else
                    button.style.backgroundColor = buttonColors[6];
            }
            if (buttonElements.has(button.innerHTML))
                buttonElements.set(button.innerHTML, button);
        }
    currentMode = newMode;
    var nFormulas = exactFormulas.length + implicitFormula.length + missingTermFormulas.length + extraTermFormulas.length;
    boldButton('1/UN', (currentMode == 'mode-unit' && unitSign == -1));
    boldButton('2nd', (currentMode == 'mode-2nd'));
    boldButton('CNST', (currentMode == 'mode-const'));
    boldButton('UNIT', (currentMode == 'mode-unit'));
    boldButton('LIST', (currentMode == 'mode-norm') && nFormulas > 0);
    boldButton('HYP', hyper);
    if (newMode == 'mode-norm' || newMode == 'mode-2nd') {
        buttonElements.get('RAD+').innerHTML = degRad;
        buttonElements.get('FLT+').innerHTML = fltSciEng;
    }
}
function fillTreeInHTML(listForTree, treeView, funcName) {
    var groups = [];
    for (var i = 0; i < listForTree.length; i++) {
        if (!groups.includes(listForTree[i].group))
            groups.push(listForTree[i].group);
    }
    var str = '';
    for (i = 0; i < groups.length; i++) {
        var symName = treeView + '-' + groups[i];
        str += `<li> <span class='caret'>${groups[i]}</span><ul id='${symName}' class='nested'> </ul> </li>`;
    }
    document.getElementById(treeView).innerHTML = str;
    for (var i = 0; i < listForTree.length; i++) {
        var unit = listForTree[i];
        var li = `<li onclick="${funcName}('${unit.name}',true)">${unit.desc} (${unit.name})</li>`;
        document.getElementById(treeView + '-' + unit.group).innerHTML += li;
    }
}
function wireUpTreeTogglerInHTML() {
    var toggler = document.getElementsByClassName('caret');
    for (var i = 0; i < toggler.length; i++) {
        toggler[i].addEventListener('click', function () {
            this.parentElement.querySelector('.nested').classList.toggle('active');
            this.classList.toggle('caret-down');
        });
    }
}
function debug() {
    var i;
    document.getElementById('operators').innerHTML = '---operators---';
    for (i = 0; i < operators.length; i++)
        document.getElementById('operators').innerHTML += '<br>' + operators[i];
    document.getElementById('operands').innerHTML = '----operands----';
    for (i = 0; i < operands.length; i++)
        document.getElementById('operands').innerHTML += '<br>' + JSON.stringify(operands[i]);
    document.getElementById('entryMode').innerHTML = '----entryMode----<BR>' + entryMode;
}
function setMessage(message) {
    if (message == undefined)
        message = 'Physical Calculator';
    document.getElementById('message').innerHTML = message;
}
function setDisplay(measure) {
    var m = measure.toDisplayForm();
    document.getElementById('value').innerHTML = m.value;
    if (m.numerator == '')
        document.getElementById('numerator').innerHTML = '';
    else
        document.getElementById('numerator').innerHTML = '&nbsp;' + m.numerator + '&nbsp;';
    if (m.denominator == '')
        document.getElementById('denominator').innerHTML = '';
    else
        document.getElementById('denominator').innerHTML = m.denominator;
}
function updateDisplay() {
    if (exponent == '')
        document.getElementById('value').innerHTML = mantisa;
    else
        document.getElementById('value').innerHTML = mantisa + '×10<sup><small>' + exponent.toString() + '</small></sup>';
}
function digitButton(symbol) {
    if (entryMode == 'number')
        operands.pop();
    if (entryMode == 'exponent') {
        if (exponent == '0')
            exponent = symbol;
        else
            exponent += symbol;
    }
    else {
        if (mantisa == '0')
            mantisa = symbol;
        else
            mantisa += symbol;
        entryMode = 'mantisa';
    }
    updateDisplay();
}
function eeButton() {
    if (exponent == '')
        exponent = '0';
    entryMode = 'exponent';
    updateDisplay();
}
function plusMinusButton() {
    if (entryMode == 'exponent') {
        if (exponent.slice(0, 1) == '-')
            exponent = exponent.slice(1);
        else {
            if (exponent.slice(0, 1) == '+')
                exponent = exponent.slice(1);
            exponent = '-' + exponent;
        }
    }
    else {
        if (mantisa.slice(0, 1) == '-')
            mantisa = mantisa.slice(1);
        else
            mantisa = '-' + mantisa;
    }
    updateDisplay();
}
function finishEntry() {
    if (entryMode == 'mantisa' || entryMode == 'exponent') {
        var disp = document.getElementById('value').innerHTML;
        disp = disp.replace('×10<sup><small>', 'E').replace('</small></sup>', '');
        operands.push(new Measurement(parseFloat(disp)));
        mantisa = '';
        exponent = '';
        entryMode = 'number';
    }
    else if (operands.length == 0)
        operands.push(new Measurement(0));
}
function precedence(op) {
    switch (op) {
        case '(':
            return 0;
        case '+':
        case '-':
            return 1;
        case '*':
        case '/':
            return 2;
        case '^':
            return 3;
    }
    return -1;
}
function evaluateTopOperator() {
    var y = operands.pop();
    var x = operands.pop();
    var ret = new Measurement(x.value, x.unitPowers, x.unitNames);
    ret.complexUnits = '';
    var i;
    var oper = operators.pop();
    switch (oper) {
        case '+':
            if (x.unitPowers.toString() != y.unitPowers.toString())
                setMessage("Can't add, units differ!");
            else
                ret.value = x.value + y.value;
            break;
        case '-':
            if (x.unitPowers.toString() != y.unitPowers.toString())
                setMessage("Can't subtract, units differ!");
            else
                ret.value = x.value - y.value;
            break;
        case '*':
            ret.value = x.value * y.value;
            for (i = 0; i < 8; i++) {
                ret.unitPowers[i] = x.unitPowers[i] + y.unitPowers[i];
                if (y.unitNames[i] != '')
                    ret.unitNames[i] = y.unitNames[i];
            }
            break;
        case '/':
            ret.value = x.value / y.value;
            for (i = 0; i < 8; i++) {
                ret.unitPowers[i] = x.unitPowers[i] - y.unitPowers[i];
                if (y.unitNames[i] != '')
                    ret.unitNames[i] = y.unitNames[i];
            }
            break;
        case 'root':
            y.value = 1 / y.value;
        case '^':
            if (y.isScalar())
                powMeasurement(ret, y.value);
            else
                setMessage('Power must be scalar!');
            break;
    }
    return ret;
}
function evaluateAtAndAbovePrecedence(prec) {
    while (operators.length > 0) {
        var op = operators[operators.length - 1];
        if (prec > precedence(op))
            break;
        else if (op == '(') {
            operators.pop();
            if (prec == 0)
                break;
        }
        else
            operands.push(evaluateTopOperator());
    }
}
function infixButton(ch) {
    finishEntry();
    switch (ch) {
        case '(':
            operators.push('(');
            break;
        case ')':
            evaluateAtAndAbovePrecedence(0);
            break;
        case '=':
            evaluateAtAndAbovePrecedence(-1);
            break;
        default:
            if (entryMode == 'infix' && operators.length > 0)
                operators.pop();
            evaluateAtAndAbovePrecedence(precedence(ch));
            entryMode = 'infix';
            operators.push(ch);
    }
    setDisplay(operands[operands.length - 1]);
}
function clearButton(all) {
    if (all) {
        operands = [];
        operands.push(new Measurement(0));
        operators = [];
        exactFormulas = [];
        implicitFormula = [];
        missingTermFormulas = [];
        extraTermFormulas = [];
        listedFormulas = implicitFormula;
        knowns = [new Measurement(0)];
        populateList();
        setMessage();
    }
    setDisplay(new Measurement(0));
    mantisa = '';
    exponent = '';
    entryMode = 'number';
    setButtonMode('mode-norm');
}
function selectUnitByName(name, tree = false) {
    finishEntry();
    var unit = findUnitByName(name);
    var top = operands[operands.length - 1];
    if (unit) {
        if (unit.isComplex() || unit.isComposite()) {
            if (!top.isScalar() && top.unitPowers.toString() != unit.units.powers.toString())
                setMessage("Units differ, can't convert!");
            else {
                if (top.isScalar()) {
                    top.value -= unit.offset;
                    top.value *= unit.factor ** unitSign;
                }
                top.unitPowers = unit.units.powers.slice();
                top.unitNames = unit.names();
            }
        }
        else {
            if (top.unitPowers[unit.index()] == 0 || (top.unitNames[unit.index()] == unit.name && top.complexUnits == ''))
                top.value *= unit.factor ** unitSign * unit.units.powers[unit.index()];
            if (top.unitNames[unit.index()] == unit.name && top.complexUnits == '')
                top.unitPowers[unit.index()] += unitSign * unit.units.powers[unit.index()];
            else if (top.unitPowers[unit.index()] == 0)
                top.unitPowers[unit.index()] = unitSign * unit.units.powers[unit.index()];
            if (top.unitPowers[unit.index()] == 0)
                top.unitNames[unit.index()] = '';
            else
                top.unitNames[unit.index()] = unit.names()[unit.index()];
        }
        top.complexUnits = (unit.isComplex() ? unit.name : '');
        setDisplay(top);
    }
    if (tree) {
        document.getElementById('unitTreeDiv').style.display = 'none';
        document.getElementById('calculator').hidden = false;
        document.getElementById('help').hidden = false;
        setButtonMode('mode-norm');
    }
}
function findConstantByName(name) {
    for (var i = 0; i < listOfConstants.length; i++)
        if (name == listOfConstants[i].name)
            return listOfConstants[i];
    return undefined;
}
function selectConstByName(button, tree = false) {
    if (entryMode == 'number')
        operands.pop();
    var cnst = findConstantByName(button);
    if (cnst) {
        operands.push(cnst.toMeasure());
        setDisplay(operands[operands.length - 1]);
    }
    if (tree) {
        document.getElementById('constTreeDiv').style.display = 'none';
        document.getElementById('calculator').hidden = false;
        document.getElementById('help').hidden = false;
        setButtonMode('mode-norm');
    }
}
function powMeasurement(top, power) {
    var undo = top.unitPowers.slice();
    for (var i = 0; i < 8; i++) {
        top.unitPowers[i] *= power;
        if (!Number.isInteger(top.unitPowers[i])) {
            top.unitPowers = undo;
            setMessage('Fractional unit powers not allowed!');
            return;
        }
    }
    top.value = top.value ** power;
}
function transcendentalOp(top, newValue) {
    if (!top.isScalarOrRadian()) {
        setMessage('Transcendental functions require scalar!');
        return;
    }
    top.unitPowers = [0, 0, 0, 0, 0, 0, 0, 0];
    top.value = newValue;
}
function unaryButton(op) {
    finishEntry();
    var top = operands[operands.length - 1];
    switch (op.toLowerCase()) {
        case '1/x':
            powMeasurement(top, -1);
            break;
        case 'x^2':
            powMeasurement(top, 2);
            break;
        case 'x^3':
            powMeasurement(top, 3);
            break;
        case 'root2':
            powMeasurement(top, 1 / 2);
            break;
        case 'root3':
            powMeasurement(top, 1 / 3);
            break;
        case 'sin':
            if (hyper)
                transcendentalOp(top, Math.sinh(top.value));
            else {
                if (degRad == 'DEG+' && top.isScalar())
                    top.value = top.value * Math.PI / 180;
                transcendentalOp(top, Math.sin(top.value));
            }
            break;
        case 'cos':
            if (hyper)
                transcendentalOp(top, Math.cosh(top.value));
            else {
                if (degRad == 'DEG+' && top.isScalar())
                    top.value = top.value * Math.PI / 180;
                transcendentalOp(top, Math.cos(top.value));
            }
            break;
        case 'tan':
            if (hyper)
                transcendentalOp(top, Math.tanh(top.value));
            else {
                if (degRad == 'DEG+' && top.isScalar())
                    top.value = top.value * Math.PI / 180;
                transcendentalOp(top, Math.tan(top.value));
            }
            break;
        case 'sin-1':
            if (hyper)
                transcendentalOp(top, Math.asinh(top.value));
            else {
                transcendentalOp(top, Math.asin(top.value));
                selectUnitByName('rad');
                if (degRad == 'DEG+')
                    selectUnitByName('deg');
            }
            break;
        case 'cos-1':
            if (hyper)
                transcendentalOp(top, Math.acosh(top.value));
            else {
                transcendentalOp(top, Math.acos(top.value));
                selectUnitByName('rad');
                if (degRad == 'DEG+')
                    selectUnitByName('deg');
            }
            break;
        case 'tan-1':
            if (hyper)
                transcendentalOp(top, Math.atanh(top.value));
            else {
                transcendentalOp(top, Math.atan(top.value));
                selectUnitByName('rad');
                if (degRad == 'DEG+')
                    selectUnitByName('deg');
            }
            break;
        case 'log':
            transcendentalOp(top, Math.log10(top.value));
            break;
        case 'ln':
            transcendentalOp(top, Math.log(top.value));
            break;
        case '10^x':
            transcendentalOp(top, 10 ** top.value);
            break;
        case 'exp':
            transcendentalOp(top, Math.exp(top.value));
            break;
    }
    setDisplay(top);
}
function findImplicitFormula() {
    var nSolutions = 0;
    var bestScore = 1E6;
    var pows = [];
    var bestSolution;
    function score() {
        var sum = 0;
        for (var i = 0; i < pows.length; i++)
            if (pows[i] == 0)
                sum += 8;
            else if (pows[i] > 0)
                sum += pows[i] * 2 - 2;
            else
                sum += pows[i] * -2 - 1;
        return sum;
    }
    function recurse(depth, sums) {
        var end = ((depth == 0) ? -1 : 4);
        for (var i = -4; i <= end; i++) {
            pows[depth] = i;
            var abs = 0;
            var sums2 = [0, 0, 0, 0, 0, 0, 0, 0];
            for (var k = 0; k < 8; k++) {
                var element = sums[k] + i * knowns[depth].unitPowers[k];
                sums2[k] = element;
                abs += Math.abs(element);
            }
            if (depth == knowns.length - 1 && abs == 0) {
                nSolutions++;
                var n = score();
                if (n < bestScore) {
                    bestSolution = pows.slice();
                    bestScore = n;
                }
            }
            else if (depth < knowns.length - 1)
                recurse(depth + 1, sums2.slice());
        }
    }
    if (knowns.length >= 3 && knowns[0].nPowers() != 0)
        recurse(0, [0, 0, 0, 0, 0, 0, 0, 0]);
    if (nSolutions == 0)
        return [];
    const roots = ['', '', 'sqrt(', 'cbrt(', 'root4('];
    var formula = 'u = ' + roots[-bestSolution[0]];
    var units = new ShortUnits(undefined, knowns[0].unitPowers);
    var vars = ['u' + '=' + units.toString()];
    var numerator = '';
    var denominator = '';
    for (var i = 1; i < knowns.length; i++) {
        var p = bestSolution[i];
        var m = 'k' + i.toString();
        var units = new ShortUnits(undefined, knowns[i].unitPowers);
        vars.push(m + '=' + units.toString());
        if (Math.abs(p) == 1)
            m += ' * ';
        else
            m += '**' + p + ' * ';
        if (p > 0)
            numerator += m;
        else if (p < 0)
            denominator += m;
    }
    numerator = numerator.substr(0, numerator.length - 3);
    denominator = denominator.substr(0, denominator.length - 3);
    if (denominator != '') {
        if (numerator == '')
            numerator = ' 1';
        formula += numerator + ' / ' + denominator;
    }
    else
        formula += numerator;
    if (bestSolution[0] < -1)
        formula += ')';
    var f = new Formulas('implicit', '', [formula], vars);
    f.matching = f.solutions[0];
    return [f];
}
function findFormula() {
    clearButton(false);
    exactFormulas = findMatchingFormulas('exact');
    missingTermFormulas = findMatchingFormulas('missing');
    extraTermFormulas = findMatchingFormulas('extra');
    implicitFormula = findImplicitFormula();
    if (exactFormulas[0])
        listedFormulas = exactFormulas;
    else
        listedFormulas = implicitFormula;
    inxFormulas = 0;
    populateList();
}
function populateList() {
    document.getElementById('pageListText').innerHTML = (inxFormulas + 1).toString() + ' of ' + listedFormulas.length.toString();
    document.getElementById('exactButton').disabled = (exactFormulas.length == 0);
    document.getElementById('exactButton').style.textDecoration = ((listedFormulas === exactFormulas) ? 'underline' : '');
    document.getElementById('missingButton').disabled = (missingTermFormulas.length == 0);
    document.getElementById('missingButton').style.textDecoration = ((listedFormulas === missingTermFormulas) ? 'underline' : '');
    document.getElementById('extraButton').disabled = (extraTermFormulas.length == 0);
    document.getElementById('extraButton').style.textDecoration = ((listedFormulas === extraTermFormulas) ? 'underline' : '');
    document.getElementById('implicitButton').disabled = (implicitFormula.length == 0);
    document.getElementById('implicitButton').style.textDecoration = ((listedFormulas === implicitFormula) ? 'underline' : '');
    document.getElementById('leftArrowButton').disabled = (inxFormulas <= 0);
    document.getElementById('rightArrowButton').disabled = (inxFormulas >= listedFormulas.length - 1);
    var formula = listedFormulas[inxFormulas];
    if (formula) {
        formula.matchVariables();
        setMessage(formula.desc + ': ' + formula.prettyMatching());
        if (listedFormulas == missingTermFormulas)
            document.getElementById('value').innerHTML = 'missing term(s)!';
        else {
            knowns[0] = new Measurement(formula.solve(), knowns[0].unitPowers, knowns[0].unitNames, knowns[0].complexUnits, knowns[0].formulaVar);
            setDisplay(knowns[0]);
        }
        document.getElementById('formula').innerHTML = formula.desc + ': ' + formula.prettyMatching();
        for (var i = 0; i < 9; i++)
            document.getElementById('list' + i.toString()).innerHTML = ((i < knowns.length) ? knowns[i].formulaVar + ' = ' + knowns[i].toDisplayForm() : '&nbsp');
    }
    else {
        var unknown = ((knowns[0].nPowers() == 0) ? '0' : '1');
        setMessage('unknown ' + unknown + ' = ' + (knowns.length - 1).toString() + ' known');
        document.getElementById('formula').innerHTML = 'no formula found';
        document.getElementById('list0').innerHTML = 'u = ' + knowns[0].toDisplayForm();
        for (var i = 1; i < 9; i++)
            document.getElementById('list' + i.toString()).innerHTML = ((i < knowns.length) ? 'k' + i.toString() + ' = ' + knowns[i].toDisplayForm() : '&nbsp');
    }
}
function toggleUnitMode() {
    finishEntry();
    if (currentMode == 'mode-unit')
        setButtonMode('mode-norm');
    else
        setButtonMode('mode-unit');
    unitSign = 1;
    finishEntry();
    setDisplay(operands[operands.length - 1]);
}
function boldButton(name, bold) {
    if (buttonElements.get(name))
        buttonElements.get(name).style.fontWeight = (bold ? 'bold' : 'normal');
}
function keyButton(evnt) {
    setMessage();
    if (!evnt)
        evnt = window.event;
    var elemt = (evnt.target || evnt.srcElement);
    var top = operands[operands.length - 1];
    switch (elemt.innerHTML.toLowerCase()) {
        case 'ce/c':
            if (entryMode == 'number' || entryMode == 'infix')
                clearButton(true);
            else
                clearButton(false);
            break;
        case 'unit':
            toggleUnitMode();
            break;
        case 'kwn':
            if (top == undefined || (top.value == 0 && top.nPowers() == 0))
                toggleUnitMode();
            else {
                knowns.push(operands.pop());
                findFormula();
            }
            break;
        case 'unkn':
            if (top == undefined || (top.value == 0 && top.nPowers() == 0))
                toggleUnitMode();
            else {
                knowns[0] = operands.pop();
                findFormula();
            }
            break;
        case 'cnst':
            setButtonMode((currentMode == 'mode-const') ? 'mode-norm' : 'mode-const');
            break;
        case 'help':
            window.open('help.html', '_blank');
            break;
        default:
            if (currentMode == 'mode-unit') {
                switch (elemt.innerHTML.toLowerCase()) {
                    case 'list':
                        document.getElementById('unitTreeDiv').style.display = 'block';
                        document.getElementById('calculator').hidden = true;
                        document.getElementById('help').hidden = true;
                        break;
                    case '1/un':
                        unitSign = -unitSign;
                        setDisplay(top);
                        break;
                    case '#':
                        top.value = top.value / top.factor();
                        top.unitPowers = [0, 0, 0, 0, 0, 0, 0, 0];
                        top.unitNames = ['', '', '', '', '', '', '', ''];
                        setDisplay(top);
                        break;
                    default:
                        selectUnitByName(elemt.innerHTML, false);
                        break;
                }
            }
            else if (currentMode == 'mode-const') {
                if (elemt.innerHTML.toLowerCase() == 'list') {
                    document.getElementById('constTreeDiv').style.display = 'block';
                    document.getElementById('calculator').hidden = true;
                    document.getElementById('help').hidden = true;
                }
                else {
                    selectConstByName(elemt.innerHTML);
                    setButtonMode('mode-norm');
                }
            }
            else {
                switch (elemt.innerHTML.toLowerCase()) {
                    case 'π':
                        selectConstByName(elemt.innerHTML);
                        break;
                    case '2nd':
                        setButtonMode((currentMode == 'mode-2nd') ? 'mode-norm' : 'mode-2nd');
                        break;
                    case 'hyp':
                        hyper = !hyper;
                        setButtonMode(currentMode);
                        break;
                    case 'deg+':
                        degRad = 'RAD+';
                        setButtonMode(currentMode);
                        if (top.isScalarOrRadian() && !top.isScalar() && top.unitNames[7] == 'deg')
                            selectUnitByName('rad');
                        break;
                    case 'rad+':
                        degRad = 'DEG+';
                        setButtonMode(currentMode);
                        if (top.isScalarOrRadian() && !top.isScalar() && top.unitNames[7] == 'rad')
                            selectUnitByName('deg');
                        break;
                    case 'flt+':
                        fltSciEng = 'SCI+';
                        setButtonMode(currentMode);
                        if (top)
                            setDisplay(top);
                        break;
                    case 'sci+':
                        fltSciEng = 'ENG+';
                        setButtonMode(currentMode);
                        if (top)
                            setDisplay(top);
                        break;
                    case 'eng+':
                        fltSciEng = 'FLT+';
                        setButtonMode(currentMode);
                        if (top)
                            setDisplay(top);
                        break;
                    case 'list':
                        document.getElementById('listDiv').style.display = 'block';
                        document.getElementById('calculator').hidden = true;
                        document.getElementById('help').hidden = true;
                        break;
                    case '0':
                    case '1':
                    case '2':
                    case '3':
                    case '4':
                    case '5':
                    case '6':
                    case '7':
                    case '8':
                    case '9':
                    case '.':
                        digitButton(elemt.innerHTML);
                        break;
                    case '+':
                    case '-':
                    case '(':
                    case ')':
                    case '=':
                        infixButton(elemt.innerHTML);
                        break;
                    case '×':
                        infixButton('*');
                        break;
                    case '÷':
                        infixButton('/');
                        break;
                    case '±':
                        plusMinusButton();
                        break;
                    case 'ee':
                        eeButton();
                        break;
                    case '¹/ₓ':
                        unaryButton('1/x');
                        break;
                    case 'ln':
                    case 'sin':
                    case 'cos':
                    case 'tan':
                    case 'log':
                        unaryButton(elemt.innerHTML);
                        break;
                    case 'eˣ':
                        unaryButton('exp');
                        break;
                    case 'yˣ':
                        infixButton('^');
                        break;
                    case '√x':
                        unaryButton('root2');
                        break;
                    case 'x²':
                        unaryButton('x^2');
                        break;
                    case 'sin⁻¹':
                        unaryButton('sin-1');
                        break;
                    case 'cos⁻¹':
                        unaryButton('cos-1');
                        break;
                    case 'tan⁻¹':
                        unaryButton('tan-1');
                        break;
                    case '10ˣ':
                        unaryButton('10^x');
                        break;
                    case 'ˣ√y':
                        infixButton('root');
                        break;
                    case '³√x':
                        unaryButton('root3');
                        break;
                    case 'x³':
                        unaryButton('x^3');
                        break;
                    case 'str':
                        finishEntry();
                        top = operands[operands.length - 1];
                        memory = top.copy();
                        break;
                    case 'rcl':
                        if (entryMode == 'number' || entryMode == 'infix') {
                            operands.push(memory);
                            setDisplay(operands[operands.length - 1]);
                        }
                        break;
                    case 'sum':
                        finishEntry();
                        top = operands[operands.length - 1];
                        if (memory.unitPowers.toString() != top.unitPowers.toString())
                            setMessage("Can't sum, units differ!");
                        else
                            memory.value += top.value;
                        break;
                    case 'exch':
                        finishEntry();
                        top = operands[operands.length - 1];
                        var temp = top.copy();
                        top = memory;
                        memory = temp;
                        setDisplay(top);
                        break;
                }
            }
    }
    if (elemt.innerHTML != '2nd' && currentMode == 'mode-2nd')
        setButtonMode('mode-norm');
    elemt.style.borderStyle = 'inset';
    elemt.style.opacity = '.5';
    var closure = function () { elemt.style.borderStyle = 'outset'; elemt.style.opacity = '1'; };
    setTimeout(closure, 200);
}
function cancelButton() {
    document.getElementById('constTreeDiv').style.display = 'none';
    document.getElementById('unitTreeDiv').style.display = 'none';
    document.getElementById('listDiv').style.display = 'none';
    document.getElementById('help').hidden = false;
    document.getElementById('calculator').hidden = false;
}
function exactListButton() {
    listedFormulas = exactFormulas;
    inxFormulas = 0;
    populateList();
}
function missingListButton() {
    listedFormulas = missingTermFormulas;
    inxFormulas = 0;
    populateList();
}
function extraListButton() {
    listedFormulas = extraTermFormulas;
    inxFormulas = 0;
    populateList();
}
function implicitListButton() {
    listedFormulas = implicitFormula;
    inxFormulas = 0;
    populateList();
}
function rightArrowListButton() {
    inxFormulas++;
    populateList();
}
function leftArrowListButton() {
    inxFormulas--;
    populateList();
}
if (navigator.userAgent.indexOf('Android') >= 0) {
    window.onscroll = function () {
        document.getElementById('page').style.height = window.innerHeight + 'px';
    };
}
function setupScroll() {
    var iphone = (navigator.userAgent.indexOf('iPhone') >= 0) || (navigator.userAgent.indexOf('iPod') >= 0);
    var ipad = (navigator.userAgent.indexOf('iPad') >= 0);
    if (iphone || ipad) {
        var height = document.documentElement.clientHeight;
        var fullscreen = (window.navigator['standalone'] == true);
        if (iphone && !fullscreen)
            height += 60;
        document.getElementById('page').style.height = height + 'px';
    }
    else if (navigator.userAgent.indexOf('Android') >= 0) {
        document.getElementById('page').style.height = (window.innerHeight + 56) + 'px';
    }
    setTimeout(scrollTo, 0, 0, 1);
}
var lastWidth = 0;
function onResize() {
    var pageWidth = document.getElementById('page').offsetWidth;
    if (lastWidth == pageWidth)
        return;
    lastWidth = pageWidth;
    setupScroll();
}
window.onresize = onResize;
document.addEventListener('keydown', handleKeydown);
function handleKeydown(e) {
    switch (e.key) {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case '.':
            digitButton(e.key);
            break;
        case '+':
        case '-':
        case '*':
        case '/':
        case '(':
        case ')':
        case '=':
            infixButton(e.key);
            break;
        case 'm':
        case 'M':
            plusMinusButton();
            break;
        case 'e':
        case 'E':
            eeButton();
            break;
        case 'Enter':
            infixButton('=');
            break;
        case 'Escape':
        case 'Esc':
            if (entryMode == 'number' || entryMode == 'infix')
                clearButton(true);
            else
                clearButton(false);
            break;
    }
}
//# sourceMappingURL=script.js.map