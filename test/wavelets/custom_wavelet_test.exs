defmodule Wavelets.CustomWaveletTest do
  use ExUnit.Case
  doctest Wavelets.CustomWavelet

  alias Wavelets.CustomWavelet

  test "creates valid custom wavelet from correct coefficients" do
    # Haar wavelet coefficients
    decomp_low = [0.7071067811865476, 0.7071067811865476]
    decomp_high = [-0.7071067811865476, 0.7071067811865476]
    recon_low = [0.7071067811865476, 0.7071067811865476]
    recon_high = [0.7071067811865476, -0.7071067811865476]

    assert {:ok, filter} =
             CustomWavelet.create(
               decomp_low,
               decomp_high,
               recon_low,
               recon_high,
               name: "custom_haar"
             )

    assert filter.name == "custom_haar"
    assert filter.support_width == 2
    assert filter.vanishing_moments == 1
  end

  test "creates valid orthogonal wavelet" do
    # Haar low-pass coefficients
    low_pass = [0.7071067811865476, 0.7071067811865476]

    assert {:ok, filter} =
             CustomWavelet.create_orthogonal(
               low_pass,
               name: "custom_orthogonal"
             )

    assert filter.name == "custom_orthogonal"
    assert filter.support_width == 2
  end

  test "creates valid biorthogonal wavelet pair" do
    # Simple biorthogonal pair
    decomp_low = [0.7071067811865476, 0.7071067811865476]
    recon_low = [0.7071067811865476, 0.7071067811865476]

    assert {:ok, filter} =
             CustomWavelet.create_biorthogonal(
               decomp_low,
               recon_low,
               name: "custom_biorthogonal"
             )

    assert filter.name == "custom_biorthogonal"
    assert length(filter.decomposition_high_pass) == length(decomp_low)
    assert length(filter.reconstruction_high_pass) == length(recon_low)
  end
end
