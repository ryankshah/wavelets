defmodule Wavelets.CustomWaveletTest do
  use ExUnit.Case
  doctest Wavelets.CustomWavelet

  alias Wavelets.CustomWavelet

  describe "create/5" do
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

    test "rejects coefficients that don't satisfy perfect reconstruction" do
      # Invalid coefficients
      coeffs = [1.0, 1.0]

      assert {:error, _} =
               CustomWavelet.create(
                 coeffs,
                 coeffs,
                 coeffs,
                 coeffs
               )
    end

    test "verifies orthogonality conditions" do
      # Non-orthogonal coefficients
      decomp_low = [1.0, 1.0]
      decomp_high = [1.0, -1.0]
      recon_low = [0.5, 0.5]
      recon_high = [0.5, -0.5]

      assert {:error, _} =
               CustomWavelet.create(
                 decomp_low,
                 decomp_high,
                 recon_low,
                 recon_high
               )
    end
  end

  describe "create_orthogonal/2" do
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
  end

  describe "create_biorthogonal/3" do
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
end
