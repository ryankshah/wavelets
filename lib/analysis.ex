defmodule Wavelets.Analysis do
  @moduledoc """
  Tools for analyzing and visualizing wavelet transforms
  """

  @doc """
  Computes energy distribution across wavelet coefficients
  """
  @spec energy_distribution(list(number) | list(list(number))) :: float | map()
  def energy_distribution(coeffs) when is_list(coeffs) do
    compute_energy(coeffs)
  end

  defp compute_energy(signal) when is_list(signal) do
    cond do
      is_number(hd(signal)) ->
        # 1D case
        signal |> Enum.map(&(&1 * &1)) |> Enum.sum()

      is_list(hd(signal)) ->
        # 2D case
        signal |> List.flatten() |> Enum.map(&(&1 * &1)) |> Enum.sum()
    end
  end
end
