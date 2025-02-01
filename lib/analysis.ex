# lib/wavelets/analysis.ex
defmodule Wavelets.Analysis do
  @moduledoc """
  Tools for analyzing and visualizing wavelet transforms
  """

  @doc """
  Computes energy distribution across wavelet coefficients
  """
  @spec energy_distribution(list(number) | list(list(number))) ::
          float | map()
  def energy_distribution(coeffs) when is_list(coeffs) do
    cond do
      is_number(hd(coeffs)) ->
        # 1D case - total signal energy
        compute_total_energy(coeffs)

      is_list(hd(coeffs)) ->
        # 2D or higher case
        compute_subband_energies(coeffs)
    end
  end

  defp compute_total_energy(signal) do
    signal
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
  end

  defp compute_subband_energies(coeffs) do
    # For 2D case, compute energies per subband
    approx_energy = compute_total_energy(List.flatten(coeffs))

    detail_energies =
      coeffs
      |> List.flatten()
      |> Enum.chunk_every(div(length(List.flatten(coeffs)), 4))
      |> Enum.map(&compute_total_energy/1)

    %{
      approximation: Enum.at(detail_energies, 0, 0.0),
      horizontal: Enum.at(detail_energies, 1, 0.0),
      vertical: Enum.at(detail_energies, 2, 0.0),
      diagonal: Enum.at(detail_energies, 3, 0.0),
      total: approx_energy
    }
  end

  @doc """
  Computes Shannon entropy of wavelet coefficients
  """
  @spec entropy(list(number) | list(list(number))) :: float
  def entropy(coeffs) when is_list(coeffs) do
    coeffs
    |> List.flatten()
    |> Enum.map(&abs/1)
    |> normalize()
    |> Enum.map(fn p ->
      if p > 0, do: -p * :math.log2(p), else: 0
    end)
    |> Enum.sum()
  end

  @doc """
  Performs statistical analysis of wavelet coefficients
  """
  @spec statistics(list(number) | list(list(number))) :: map()
  def statistics(coeffs) when is_list(coeffs) do
    flat_coeffs = List.flatten(coeffs)

    %{
      mean: mean(flat_coeffs),
      variance: variance(flat_coeffs),
      skewness: skewness(flat_coeffs),
      kurtosis: kurtosis(flat_coeffs),
      max_amplitude: Enum.max_by(flat_coeffs, &abs/1),
      sparsity: sparsity(flat_coeffs)
    }
  end

  # Private helper functions
  defp compute_energy_1d(coeffs) do
    coeffs
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
  end

  defp compute_energy_2d(coeffs) do
    total_energy = compute_energy_1d(List.flatten(coeffs))
    level_count = trunc(:math.log2(length(coeffs)))

    subbands =
      for level <- 1..level_count do
        size = trunc(:math.pow(2, level_count - level))

        {
          approximation_energy(coeffs, size),
          horizontal_energy(coeffs, size),
          vertical_energy(coeffs, size),
          diagonal_energy(coeffs, size)
        }
      end

    %{
      total_energy: total_energy,
      subbands: subbands
    }
  end

  defp approximation_energy(coeffs, size) do
    coeffs
    |> Enum.take(size)
    |> Enum.map(&Enum.take(&1, size))
    |> List.flatten()
    |> compute_energy_1d()
  end

  defp horizontal_energy(coeffs, size) do
    coeffs
    |> Enum.take(size)
    |> Enum.map(&Enum.slice(&1, size..-1))
    |> List.flatten()
    |> compute_energy_1d()
  end

  defp vertical_energy(coeffs, size) do
    coeffs
    |> Enum.slice(size..-1)
    |> Enum.map(&Enum.take(&1, size))
    |> List.flatten()
    |> compute_energy_1d()
  end

  defp diagonal_energy(coeffs, size) do
    coeffs
    |> Enum.slice(size..-1)
    |> Enum.map(&Enum.slice(&1, size..-1))
    |> List.flatten()
    |> compute_energy_1d()
  end

  defp normalize(values) do
    sum = Enum.sum(Enum.map(values, &abs/1))

    if sum > 0 do
      Enum.map(values, &(abs(&1) / sum))
    else
      values
    end
  end

  defp mean(values) do
    Enum.sum(values) / length(values)
  end

  defp variance(values) do
    m = mean(values)

    values
    |> Enum.map(&:math.pow(&1 - m, 2))
    |> mean()
  end

  defp skewness(values) do
    m = mean(values)
    v = variance(values)
    std = :math.sqrt(v)

    if std > 0 do
      values
      |> Enum.map(&:math.pow((&1 - m) / std, 3))
      |> mean()
    else
      0.0
    end
  end

  defp kurtosis(values) do
    m = mean(values)
    v = variance(values)
    std = :math.sqrt(v)

    if std > 0 do
      values
      |> Enum.map(&:math.pow((&1 - m) / std, 4))
      |> mean()
      # Excess kurtosis
      |> Kernel.-(3)
    else
      0.0
    end
  end

  defp sparsity(values) do
    max_amp = Enum.max_by(values, &abs/1)
    threshold = 0.001 * max_amp
    count_small = Enum.count(values, &(abs(&1) < threshold))
    count_small / length(values)
  end

  @doc """
  Generates visualization data for plotting
  """
  @spec visualization_data(list(number) | list(list(number))) :: map()
  def visualization_data(coeffs) when is_list(coeffs) do
    cond do
      is_number(hd(coeffs)) ->
        # 1D visualization data
        generate_1d_viz_data(coeffs)

      is_list(hd(coeffs)) ->
        # 2D visualization data
        generate_2d_viz_data(coeffs)
    end
  end

  defp generate_1d_viz_data(coeffs) do
    %{
      values: coeffs,
      x_range: 0..(length(coeffs) - 1) |> Enum.to_list(),
      type: "1d",
      energy: compute_energy_1d(coeffs)
    }
  end

  defp generate_2d_viz_data(coeffs) do
    %{
      values: coeffs,
      dimensions: {length(hd(coeffs)), length(coeffs)},
      type: "2d",
      energy: compute_energy_2d(coeffs)
    }
  end

  @doc """
  Computes relative energy distribution across levels
  """
  @spec level_energy_distribution(list(number) | list(list(number))) ::
          list({integer, float})
  def level_energy_distribution(coeffs) when is_list(coeffs) do
    total_energy = compute_energy_1d(List.flatten(coeffs))

    if total_energy > 0 do
      coeffs_by_level = split_by_level(coeffs)

      coeffs_by_level
      |> Enum.with_index()
      |> Enum.map(fn {level_coeffs, level} ->
        energy = compute_energy_1d(List.flatten(level_coeffs))
        {level, energy / total_energy}
      end)
    else
      []
    end
  end

  defp split_by_level(coeffs) when is_list(coeffs) do
    cond do
      is_number(hd(coeffs)) ->
        # Base case: single level
        [coeffs]

      is_list(hd(coeffs)) ->
        # Split into approximation and details
        n = length(coeffs)
        half_n = div(n, 2)

        approx = Enum.take(coeffs, half_n)
        details = Enum.drop(coeffs, half_n)

        # Recursively split approximation
        split_by_level(approx) ++ [details]
    end
  end
end
