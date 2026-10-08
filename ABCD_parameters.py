"""
ABCD Parameters Calculator
--------------------------
A menu-driven program for engineering students to calculate the ABCD
(transmission) parameters of a two-port network such as a transmission line,
to cascade networks, and to find sending-end and receiving-end quantities.

Two-port equations:
    Vs = A * Vr + B * Ir
    Is = C * Vr + D * Ir

    | Vs |   | A  B |   | Vr |
    |    | = |      | * |    |
    | Is |   | C  D |   | Ir |

Basic elements:
    Series impedance Z      [ 1  Z ]          Shunt admittance Y     [ 1  0 ]
                            [ 0  1 ]                                 [ Y  1 ]
    Ideal transformer (a = N1/N2)   [ a   0  ]
                                    [ 0  1/a ]

Transmission line models (Z = total series impedance, Y = total shunt admittance):
    Short line     A = D = 1          B = Z                C = 0
    Nominal-pi     A = D = 1 + ZY/2   B = Z                C = Y(1 + ZY/4)
    Nominal-T      A = D = 1 + ZY/2   B = Z(1 + ZY/4)     C = Y
    Long line      A = D = cosh(gl)   B = Zc*sinh(gl)      C = sinh(gl)/Zc

Properties:
    Reciprocal network:   A*D - B*C = 1
    Symmetrical network:  A = D
    Cascade of networks:  multiply the ABCD matrices in order (source to load)

Inverse relations (receiving end from sending end):
    Vr = D * Vs - B * Is
    Ir = -C * Vs + A * Is

Conversion to other parameters:
    Z-parameters:  Z11 = A/C,  Z12 = (AD - BC)/C,  Z21 = 1/C,  Z22 = D/C
    Y-parameters:  Y11 = D/B,  Y12 = -(AD - BC)/B, Y21 = -1/B, Y22 = A/B
"""

import cmath
import math

SQRT3 = math.sqrt(3)

MODELS = {
    "1": ("short", "Short line"),
    "2": ("pi", "Medium line (nominal-pi)"),
    "3": ("t", "Medium line (nominal-T)"),
    "4": ("long", "Long line (distributed)"),
}


# ---------------------------------------------------------------- input helpers
def get_positive_float(prompt):
    """Keep asking until the user enters a valid positive number."""
    while True:
        try:
            value = float(input(prompt))
            if value <= 0:
                print("  Please enter a value greater than zero.")
                continue
            return value
        except ValueError:
            print("  Invalid input. Please enter a number.")


def get_nonneg_float(prompt):
    """Ask for a number that may be zero but not negative."""
    while True:
        try:
            value = float(input(prompt))
            if value < 0:
                print("  Please enter zero or a positive value.")
                continue
            return value
        except ValueError:
            print("  Invalid input. Please enter a number.")


def get_float(prompt):
    """Ask for any number (zero and negative values allowed)."""
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("  Invalid input. Please enter a number.")


def get_power_factor(prompt="Power factor (0 to 1): "):
    """Ask for a power factor between 0 (exclusive) and 1 (inclusive)."""
    while True:
        pf = get_positive_float(prompt)
        if pf <= 1:
            return pf
        print("  Power factor cannot be greater than 1.")


def get_leading():
    """Ask whether the power factor is lagging or leading."""
    while True:
        t = input("Power factor type - (L)agging or (E) leading: ").strip().lower()
        if t in ("l", "lagging"):
            return False
        if t in ("e", "leading"):
            return True
        print("  Please enter L or E.")


def ask_model():
    """Ask which transmission line model to use."""
    for key, (_, name) in MODELS.items():
        print(f"   {key}. {name}")
    while True:
        pick = input("  Select model (1-4): ").strip()
        if pick in MODELS:
            return MODELS[pick][0]
        print("  Please enter a number from 1 to 4.")


def get_line_data(model):
    """Ask for the line parameters (per phase, per km)."""
    length = get_positive_float("Line length (km): ")
    r = get_positive_float("Resistance r (ohm/km): ")
    x = get_positive_float("Inductive reactance x (ohm/km): ")
    if model == "long":
        b = get_positive_float("Shunt susceptance b (microsiemens/km): ")
    else:
        b = get_nonneg_float("Shunt susceptance b (microsiemens/km, 0 if ignored): ")
    return length, r, x, b


# ------------------------------------------------------------------ core maths
def series_impedance(z):
    """ABCD of a series impedance."""
    return 1 + 0j, z, 0j, 1 + 0j


def shunt_admittance(y):
    """ABCD of a shunt admittance."""
    return 1 + 0j, 0j, y, 1 + 0j


def ideal_transformer(a):
    """ABCD of an ideal transformer with turns ratio a = N1/N2."""
    return complex(a), 0j, 0j, complex(1 / a)


def multiply(m1, m2):
    """Cascade two networks: m1 is nearer the source, m2 nearer the load."""
    a1, b1, c1, d1 = m1
    a2, b2, c2, d2 = m2
    return (a1 * a2 + b1 * c2,
            a1 * b2 + b1 * d2,
            c1 * a2 + d1 * c2,
            c1 * b2 + d1 * d2)


def determinant(m):
    """AD - BC (equals 1 for a reciprocal network)."""
    a, b, c, d = m
    return a * d - b * c


def line_abcd(model, length, r, x, b_us):
    """Return the ABCD constants of a transmission line model."""
    z_per_km = complex(r, x)
    y_per_km = complex(0, b_us * 1e-6)
    z_total = z_per_km * length
    y_total = y_per_km * length

    if model == "short":
        return 1 + 0j, z_total, 0j, 1 + 0j
    if model == "pi":
        a = 1 + z_total * y_total / 2
        return a, z_total, y_total * (1 + z_total * y_total / 4), a
    if model == "t":
        a = 1 + z_total * y_total / 2
        return a, z_total * (1 + z_total * y_total / 4), y_total, a

    gamma = cmath.sqrt(z_per_km * y_per_km)
    zc = cmath.sqrt(z_per_km / y_per_km)
    gl = gamma * length
    a = cmath.cosh(gl)
    return a, zc * cmath.sinh(gl), cmath.sinh(gl) / zc, a


# ------------------------------------------------------------------- display
def fmt_polar(value, digits=4):
    mag, ang = cmath.polar(value)
    return f"{mag:.{digits}f} < {math.degrees(ang):.2f} deg"


def show_abcd(m):
    """Print the ABCD constants and the standard checks."""
    a, b, c, d = m
    print(f"  A = {fmt_polar(a)}")
    print(f"  B = {fmt_polar(b)} ohm")
    print(f"  C = {fmt_polar(c, 6)} S")
    print(f"  D = {fmt_polar(d)}")
    det = determinant(m)
    print(f"\n  Check A*D - B*C = {det.real:.5f} {'+' if det.imag >= 0 else '-'} j{abs(det.imag):.5f}  (should be 1)")
    print(f"  Symmetrical network (A = D): {'Yes' if cmath.isclose(a, d, rel_tol=1e-6) else 'No'}")


def line_abcd_option():
    print("  Choose the line model:")
    model = ask_model()
    length, r, x, b = get_line_data(model)
    m = line_abcd(model, length, r, x, b)

    print("\n  ----- ABCD Parameters -----")
    show_abcd(m)

    # Z and Y parameters
    a, b_, c, d = m
    print("\n  ----- Equivalent Z and Y parameters -----")
    if abs(c) < 1e-15:
        print("  C = 0, so the Z-parameters do not exist for this model.")
    else:
        print(f"  Z11 = {fmt_polar(a / c, 2)} ohm    Z12 = {fmt_polar(determinant(m) / c, 2)} ohm")
        print(f"  Z21 = {fmt_polar(1 / c, 2)} ohm    Z22 = {fmt_polar(d / c, 2)} ohm")
    print(f"  Y11 = {fmt_polar(d / b_, 6)} S    Y12 = {fmt_polar(-determinant(m) / b_, 6)} S")
    print(f"  Y21 = {fmt_polar(-1 / b_, 6)} S    Y22 = {fmt_polar(a / b_, 6)} S")


def cascade_option():
    print("  Build the network from the SOURCE end to the LOAD end.")
    total = None
    count = 0
    while True:
        print("\n  Add an element:")
        print("   1. Series impedance")
        print("   2. Shunt admittance")
        print("   3. Ideal transformer")
        print("   4. Transmission line section")
        print("   0. Finish")
        pick = input("  Enter your choice: ").strip()

        if pick == "1":
            r = get_nonneg_float("  Resistance R (ohm): ")
            x = get_float("  Reactance X (ohm, - for capacitive): ")
            element = series_impedance(complex(r, x))
        elif pick == "2":
            g = get_nonneg_float("  Conductance G (mS): ")
            b = get_float("  Susceptance B (mS, + capacitive, - inductive): ")
            element = shunt_admittance(complex(g, b) * 1e-3)
        elif pick == "3":
            element = ideal_transformer(get_positive_float("  Turns ratio a = N1/N2: "))
        elif pick == "4":
            print("  Choose the line model:")
            model = ask_model()
            length, r, x, b = get_line_data(model)
            element = line_abcd(model, length, r, x, b)
        elif pick == "0":
            if total is None:
                print("  Add at least one element first.")
                continue
            break
        else:
            print("  Invalid choice.")
            continue

        total = element if total is None else multiply(total, element)
        count += 1
        print(f"  Element {count} added.")

    print(f"\n  ----- Overall ABCD of {count} cascaded element(s) -----")
    show_abcd(total)


def sending_end_option():
    print("  Choose the line model:")
    model = ask_model()
    length, r, x, b = get_line_data(model)
    v_kv = get_positive_float("Receiving-end line voltage (kV): ")
    p_mw = get_positive_float("Receiving-end load power (MW, three-phase): ")
    pf = get_power_factor("Load power factor (0 to 1): ")
    leading = get_leading() if pf < 1 else False

    a, b_, c, d = line_abcd(model, length, r, x, b)
    vr = complex(v_kv * 1000 / SQRT3, 0)
    ir_mag = p_mw * 1e6 / (SQRT3 * v_kv * 1000 * pf)
    ir = cmath.rect(ir_mag, math.acos(pf) * (1 if leading else -1))

    vs = a * vr + b_ * ir
    i_s = c * vr + d * ir
    regulation = (abs(vs) / abs(a) - abs(vr)) / abs(vr) * 100
    p_send = 3 * (vs * i_s.conjugate()).real

    print("\n  ----- Sending End from Receiving End -----")
    print(f"  Receiving-end current   = {ir_mag:.2f} A")
    print(f"  Sending-end line voltage = {SQRT3 * abs(vs) / 1000:.3f} kV (angle {math.degrees(cmath.phase(vs)):.2f} deg)")
    print(f"  Sending-end current     = {abs(i_s):.2f} A")
    print(f"  Sending-end power factor = {math.cos(cmath.phase(vs) - cmath.phase(i_s)):.4f}")
    print(f"  Sending-end power       = {p_send / 1e6:.3f} MW")
    print(f"  Voltage regulation      = {regulation:.3f} %")
    print(f"  Efficiency              = {p_mw * 1e6 / p_send * 100:.2f} %")


def receiving_end_option():
    print("  Choose the line model:")
    model = ask_model()
    length, r, x, b = get_line_data(model)
    vs_kv = get_positive_float("Sending-end line voltage (kV): ")
    is_amp = get_positive_float("Sending-end current (A): ")
    pf = get_power_factor("Sending-end power factor (0 to 1): ")
    leading = get_leading() if pf < 1 else False

    a, b_, c, d = line_abcd(model, length, r, x, b)
    vs = complex(vs_kv * 1000 / SQRT3, 0)
    i_s = cmath.rect(is_amp, math.acos(pf) * (1 if leading else -1))

    vr = d * vs - b_ * i_s
    ir = -c * vs + a * i_s
    regulation = (abs(vs) / abs(a) - abs(vr)) / abs(vr) * 100
    p_recv = 3 * (vr * ir.conjugate()).real
    p_send = 3 * (vs * i_s.conjugate()).real

    print("\n  ----- Receiving End from Sending End -----")
    print(f"  Receiving-end line voltage = {SQRT3 * abs(vr) / 1000:.3f} kV (angle {math.degrees(cmath.phase(vr)):.2f} deg)")
    print(f"  Receiving-end current     = {abs(ir):.2f} A")
    print(f"  Receiving-end power factor = {math.cos(cmath.phase(vr) - cmath.phase(ir)):.4f}")
    print(f"  Receiving-end power       = {p_recv / 1e6:.3f} MW")
    print(f"  Voltage regulation        = {regulation:.3f} %")
    print(f"  Efficiency                = {p_recv / p_send * 100:.2f} %")


def menu():
    print("\n" + "=" * 56)
    print("          ABCD PARAMETERS CALCULATOR")
    print("=" * 56)
    print(" 1. ABCD parameters of a transmission line (+ Z, Y)")
    print(" 2. Cascade of networks (series, shunt, transformer, line)")
    print(" 3. Sending end from receiving end")
    print(" 4. Receiving end from sending end (inverse)")
    print(" 0. Exit")
    print("-" * 56)


def main():
    while True:
        menu()
        choice = input("Enter your choice: ").strip()

        if choice == "1":
            line_abcd_option()
        elif choice == "2":
            cascade_option()
        elif choice == "3":
            sending_end_option()
        elif choice == "4":
            receiving_end_option()
        elif choice == "0":
            print("\nThank you for using the calculator. Goodbye!")
            break
        else:
            print("  Invalid choice. Please select from the menu.")

