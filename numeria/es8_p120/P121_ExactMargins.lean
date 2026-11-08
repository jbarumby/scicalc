import Mathlib.Tactic
/-! NUMERIA P121: original P111/P115 margin arithmetic on exact rational witnesses.
    All six theorems below verify rational identities only.
    They do not independently prove Bareiss/FLINT or Arb analytic inclusion. -/
namespace Numeria.ES8P121
def shift_p111 : ℚ := (7354363380645 : ℚ) / 154742504910672534362390528
def rho_p111 : ℚ := (16745517027 : ℚ) / 77371252455336267181195264
def lower_p111 : ℚ := (7320872346591 : ℚ) / 154742504910672534362390528
theorem identity_p111 : lower_p111 = shift_p111 - rho_p111 := by
  norm_num [lower_p111, shift_p111, rho_p111]
theorem margin_p111 : rho_p111 < shift_p111 := by
  norm_num [rho_p111, shift_p111]
theorem lower_positive_p111 : 0 < lower_p111 := by
  norm_num [lower_p111]
def shift_p115 : ℚ := (27215954001 : ℚ) / 604462909807314587353088
def rho_p115 : ℚ := (261648779 : ℚ) / 1208925819614629174706176
def lower_p115 : ℚ := (54170259223 : ℚ) / 1208925819614629174706176
theorem identity_p115 : lower_p115 = shift_p115 - rho_p115 := by
  norm_num [lower_p115, shift_p115, rho_p115]
theorem margin_p115 : rho_p115 < shift_p115 := by
  norm_num [rho_p115, shift_p115]
theorem lower_positive_p115 : 0 < lower_p115 := by
  norm_num [lower_p115]
end Numeria.ES8P121
