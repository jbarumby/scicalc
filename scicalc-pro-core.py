# SciCalc Pro — motor simbólico isolado (Pyodide/SymPy)
# © 2026 João Carlos Barumby. Integração adicional ao SciCalc.
# O navegador executa este arquivo num Web Worker; nenhum cálculo vai a um servidor.
import json
import re
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, convert_xor, implicit_multiplication_application
from sympy.calculus.util import continuous_domain

_X = sp.Symbol('x', real=True)
_FUNCS = {
    'x': _X, 'pi': sp.pi, 'E': sp.E, 'e': sp.E,
    'oo': sp.oo, 'inf': sp.oo,
    'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan,
    'asin': sp.asin, 'acos': sp.acos, 'atan': sp.atan,
    'sinh': sp.sinh, 'cosh': sp.cosh, 'tanh': sp.tanh,
    'exp': sp.exp, 'log': sp.log, 'ln': sp.log, 'sqrt': sp.sqrt,
    'abs': sp.Abs, 'Abs': sp.Abs, 'erf': sp.erf,
    'gamma': sp.gamma,
}
_GLOBALS = {
    '__builtins__': {}, 'Integer': sp.Integer, 'Float': sp.Float,
    'Rational': sp.Rational, 'Symbol': sp.Symbol,
    'Add': sp.Add, 'Mul': sp.Mul, 'Pow': sp.Pow,
}
_TRANSFORMS = standard_transformations + (convert_xor, implicit_multiplication_application)
_IDENTIFIERS = set(_FUNCS)

def _parse(value, label='expressão'):
    source = str(value).strip()
    if not source or len(source) > 240:
        raise ValueError(f'{label}: digite até 240 caracteres.')
    if not re.fullmatch(r'[A-Za-z0-9+\-*/^()., \t]+', source):
        raise ValueError(f'{label}: use apenas números, x e funções matemáticas permitidas.')
    tokens = re.findall(r'[A-Za-z][A-Za-z0-9]*', source)
    if any(token not in _IDENTIFIERS for token in tokens):
        unknown = sorted(set(token for token in tokens if token not in _IDENTIFIERS))
        raise ValueError('Identificador não permitido: ' + ', '.join(unknown))
    if re.search(r'(\*\*|\^)\s*[0-9]{4,}', source):
        raise ValueError('Expoente excessivamente grande.')
    try:
        result = parse_expr(source, local_dict=_FUNCS, global_dict=_GLOBALS,
                            transformations=_TRANSFORMS, evaluate=True)
    except Exception as exc:
        raise ValueError(f'Não foi possível interpretar {label}: {type(exc).__name__}.')
    if not isinstance(result, sp.Basic):
        raise ValueError('A entrada não representa uma expressão matemática.')
    if result.free_symbols - {_X}:
        raise ValueError('Esta versão aceita apenas a variável x.')
    return result

def _tex(value):
    return sp.latex(value)

def _numeric(value):
    if value in (sp.oo, -sp.oo, sp.zoo, sp.nan):
        return None
    try:
        estimate = sp.N(value, 16)
        if estimate.is_real and estimate.is_finite:
            return str(estimate)
    except Exception:
        pass
    return None

def _result(expr, answer, status, note, **extra):
    output = {
        'ok': True, 'input_tex': _tex(expr),
        'result_tex': _tex(answer), 'result_text': str(answer),
        'approx': _numeric(answer), 'status': status, 'note': note
    }
    output.update(extra)
    return output

def scicalc_pro_solve(request_json):
    req = json.loads(request_json)
    action = str(req.get('action', ''))
    expr = _parse(req.get('expr', ''), 'função')
    if action == 'derivative':
        order = int(req.get('order', 1))
        if not 1 <= order <= 6:
            raise ValueError('Escolha uma ordem de derivação entre 1 e 6.')
        value = sp.diff(expr, _X, order)
        output = _result(expr, value, 'Simbólico', f'Derivada de ordem {order}, calculada por diferenciação simbólica.')
        point = str(req.get('point', '')).strip()
        if point:
            p = _parse(point, 'ponto de avaliação')
            if p.free_symbols or p.is_real is not True:
                raise ValueError('O ponto deve ser um valor real.')
            evaluated = value.subs(_X, p)
            output['point_tex'] = _tex(p)
            output['point_result_tex'] = _tex(evaluated)
            output['point_approx'] = _numeric(evaluated)
        return json.dumps(output, ensure_ascii=False)

    if action == 'integral':
        low = str(req.get('lower', '')).strip()
        high = str(req.get('upper', '')).strip()
        if bool(low) != bool(high):
            raise ValueError('Informe ambos os limites ou deixe os dois em branco.')
        if not low:
            answer = sp.integrate(expr, _X)
            if answer.has(sp.Integral):
                output = _result(expr, answer, 'Não resolvida em forma fechada',
                  'O integrador não encontrou uma primitiva fechada. A integral formal permanece indicada.')
            else:
                output = _result(expr, answer, 'Simbólico',
                  'Primitiva simbólica. Some a constante de integração C.')
                # Uma identidade algébrica nula permite indicar verificação, sem inferir além dela.
                try:
                    output['verified'] = bool(sp.simplify(sp.diff(answer, _X) - expr) == 0)
                except Exception:
                    output['verified'] = False
            output['constant'] = not answer.has(sp.Integral)
            return json.dumps(output, ensure_ascii=False)

        a, b = _parse(low, 'limite inferior'), _parse(high, 'limite superior')
        if a.free_symbols or b.free_symbols:
            raise ValueError('Os limites devem ser números ou ±inf.')
        if (a not in (sp.oo, -sp.oo) and a.is_real is False) or (b not in (sp.oo, -sp.oo) and b.is_real is False):
            raise ValueError('Use limites reais.')
        integral = sp.Integral(expr, (_X, a, b))
        answer = sp.integrate(expr, (_X, a, b))
        if answer.has(sp.Integral):
            output = _result(integral, answer, 'Integral formal',
                  'Sem avaliação exata conclusiva. Evitamos apresentar uma aproximação como prova de convergência.')
            # Aproximação apenas em intervalo compacto onde a continuidade foi confirmada.
            if a.is_finite and b.is_finite and a.is_number and b.is_number:
                try:
                    interval = sp.Interval(sp.Min(a, b), sp.Max(a, b))
                    domain = continuous_domain(expr, _X, interval)
                    if domain == interval:
                        numeric = sp.N(integral, 15)
                        if numeric.is_real and numeric.is_finite:
                            output['approx'] = str(numeric)
                            output['note'] = 'Estimativa numérica em intervalo de continuidade; não é certificação rigorosa do erro.'
                except Exception:
                    pass
            return json.dumps(output, ensure_ascii=False)
        return json.dumps(_result(integral, answer, 'Exato / simbólico',
            'Resultado simbólico quando resolúvel; divergências e singularidades devem ser interpretadas no domínio real.'), ensure_ascii=False)

    if action == 'limit':
        target = _parse(req.get('target', '0'), 'ponto do limite')
        if target.free_symbols or (target not in (sp.oo, -sp.oo) and target.is_real is False):
            raise ValueError('O ponto do limite precisa ser real ou ±inf.')
        direction = str(req.get('direction', 'both'))
        if direction not in ('both', 'left', 'right'):
            raise ValueError('Direção inválida.')
        if target in (sp.oo, -sp.oo):
            direction = 'right'
        if direction == 'both':
            left = sp.limit(expr, _X, target, dir='-')
            right = sp.limit(expr, _X, target, dir='+')
            same = (left == right)
            if not same:
                try:
                    same = sp.simplify(left - right) == 0
                except Exception:
                    pass
            if not same:
                output = _result(expr, sp.Symbol('DNE'), 'Não existe (bilateral)',
                    'Os limites laterais são diferentes. Não há limite bilateral.',
                    left_tex=_tex(left), right_tex=_tex(right))
                output['approx'] = None
                return json.dumps(output, ensure_ascii=False)
            answer = left
        else:
            answer = sp.limit(expr, _X, target, dir='-' if direction == 'left' else '+')
        if isinstance(answer, sp.AccumBounds) or answer.has(sp.Limit) or answer in (sp.nan, sp.zoo):
            status, note = 'Inconclusivo / não convergente', 'Não foi encontrado um limite pontual finito ou infinito bem definido.'
        else:
            status, note = 'Limite simbólico', 'O resultado usa a análise simbólica de limites, sem estimativa por amostragem.'
        output = _result(expr, answer, status, note)
        output['target_tex'] = _tex(target)
        output['direction'] = direction
        return json.dumps(output, ensure_ascii=False)
    raise ValueError('Operação não reconhecida.')
